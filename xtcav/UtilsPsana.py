import numpy as np
import psana
import warnings
import time
from Utils import ROIMetrics, GlobalCalibration, ShotToShotParameters
import Constants

def GetCameraSaturationValue(evt):
    try:
        analysis_version = psana.Detector(Constants.ANALYSIS_VERSION)
        if analysis_version(evt) is not None:
            return (1<<12)-1
    except:
        pass

    return (1<<14)-1
    

def GetGlobalXTCAVCalibration(evt):
    """
    Obtain the global XTCAV calibration form the epicsStore
    Arguments:
      epicsStore
    Output:
      globalCalibration: struct with the parameters
      ok: if all the data was retrieved correctly
    """
    umperpix = psana.Detector('XTCAV_calib_umPerPx')
    strstrength = psana.Detector('XTCAV_strength_par_S')
    rfampcalib = psana.Detector('XTCAV_Amp_Des_calib_MV')
    rfphasecalib = psana.Detector('XTCAV_Phas_Des_calib_deg')
    dumpe = psana.Detector('XTCAV_Beam_energy_dump_GeV')
    dumpdisp = psana.Detector('XTCAV_calib_disp_posToEnergy')

    global_calibration = GlobalCalibration(
        umperpix=umperpix(evt), 
        strstrength=strstrength(evt), 
        rfampcalib=rfampcalib(evt), 
        rfphasecalib=rfphasecalib(evt), 
        dumpe=dumpe(evt), 
        dumpdisp=dumpdisp(evt)
    )
        
    for k,v in global_calibration._asdict().iteritems():
        if not v:
            warnings.warn_explicit('No XTCAV Calibration for epics variable ' + k, UserWarning,'XTCAV',0)
            return None

    return global_calibration
                          

def GetXTCAVImageROI(evt):

    roiXN=psana.Detector(Constants.ROI_SIZE_X)
    roiX=psana.Detector(Constants.ROI_START_X)
    roiYN=psana.Detector(Constants.ROI_SIZE_Y)
    roiY=psana.Detector(Constants.ROI_START_Y)

    xN = roiXN(evt)  #Size of the image in X                           
    x0 = roiX(evt)    #Position of the first pixel in x
    yN = roiYN(evt)  #Size of the image in Y 
    y0 = roiY(evt)    #Position of the first pixel in y
    
    if not xN:       
        warnings.warn_explicit('No XTCAV ROI info',UserWarning,'XTCAV',0)
        return None
        
    x = x0+np.arange(0, xN) 
    y = y0+np.arange(0, yN) 

    return ROIMetrics(xN, x0, yN, y0, x, y) 


def GetShotToShotParameters(ebeam, gasdetector, evt_id):
    time = evt_id.time()
    sec  = time[0]
    nsec = time[1]
    unixtime = int((sec<<32)|nsec)
    fiducial = evt_id.fiducials()

    energydetector = Constants.ENERGY_DETECTOR
 
    if ebeam:    
        ebeamcharge=ebeam.ebeamCharge()
        xtcavrfamp=ebeam.ebeamXTCAVAmpl()
        xtcavrfphase=ebeam.ebeamXTCAVPhase()
        dumpecharge=ebeam.ebeamDumpCharge()*Constants.E_CHARGE #In C 
        
        if gasdetector:
            energydetector=(gasdetector.f_11_ENRC()+gasdetector.f_12_ENRC())/2 
            return ShotToShotParameters(ebeamcharge = ebeamcharge, 
                xtcavrfphase = xtcavrfphase, xtcavrfamp = xtcavrfamp, 
                dumpecharge = dumpecharge, xrayenergy = 1e-3*energydetector, 
                unixtime = unixtime, fiducial = fiducial)     
        else:   #Some hardcoded values
            warnings.warn_explicit('No gas detector info',UserWarning,'XTCAV',0)
                
    else:    
        warnings.warn_explicit('No ebeamv info',UserWarning,'XTCAV',0)
    
    return ShotToShotParameters(unixtime = unixtime, fiducial = fiducial, valid = 0)
        


def DivideImageTasks(first_image, last_image, rank, size):
    num_shots = last_image - first_image
    if num_shots <= 0:
        return np.empty()
    tiling = np.arange(rank*4, rank*4+4,1) #  returns [0, 1, 2, 3] if e.g. rank == 0 and size == 4:
    comb1 = np.tile(tiling, np.ceil(num_shots/(4.*size)).astype(int))  # returns [0, 1, 2, 3, 0, 1, 2, 3, ...]        
    comb2 = np.repeat(np.arange(0, np.ceil(num_shots/(4.*size)), 1), 4) # returns [0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, ...]
            #  list of shot numbers assigned to this core
    main = comb2*4*size + comb1  + first_image # returns [  0.   1.   2.   3.  16.  17.  18.  19.  32.  33. ... ]
    main = np.delete(main, np.where(main>=last_image) )  # remove element if greater or equal to maximum number of shots in run
    return main.astype(int)

