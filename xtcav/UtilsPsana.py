import numpy as np
import psana
import warnings
import time
from Metrics import *
import Constants

def GetCameraSaturationValue(run, xtcav_camera, start=None):
    analysis_version = psana.Detector('XTCAV_Analysis_Version')
    times = run.times()

    end_of_images = len(times)
    if not start in range(end_of_images + 1):
        start = 0

    for t in range(start, end_of_images):
        evt=run.event(times[t])
        img = xtcav_camera.image(evt)
        # skip if empty image
        if img is None: 
            continue

        if analysis_version(evt) is not None:
            return (1<<12)-1, t
        else:
            return (1<<14)-1, t

    return None, end_of_images
    

def GetGlobalXTCAVCalibration(run, xtcav_camera, start=None):
    """
    Obtain the global XTCAV calibration form the epicsStore
    Arguments:
      epicsStore
    Output:
      globalCalibration: struct with the parameters
      ok: if all the data was retrieved correctly
    """

    umperpix=psana.Detector('XTCAV_calib_umPerPx')
    strstrength=psana.Detector('XTCAV_strength_par_S')
    rfampcalib=psana.Detector('XTCAV_Amp_Des_calib_MV')
    rfphasecalib=psana.Detector('XTCAV_Phas_Des_calib_deg')
    dumpe=psana.Detector('XTCAV_Beam_energy_dump_GeV')
    dumpdisp=psana.Detector('XTCAV_calib_disp_posToEnergy')
    times = run.times()
    
    end_of_images = len(times)
    if not start in range(end_of_images + 1):
        start = 0

    for t in range(start, end_of_images):
        evt = run.event(times[t])
        img = xtcav_camera.image(evt)
        # skip if empty image
        if img is None: 
            continue

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
                continue

        return global_calibration, t
                
    return None, end_of_images
          

def GetXTCAVImageROI(run, xtcav_camera, start=None):

    roiXN=psana.Detector('XTCAV_ROI_sizeX')
    roiX=psana.Detector('XTCAV_ROI_startX')
    roiYN=psana.Detector('XTCAV_ROI_sizeY')
    roiY=psana.Detector('XTCAV_ROI_startY')
    times = run.times()

    end_of_images = len(times)
    if not start in range(end_of_images + 1):
        start = 0

    for t in range(start, end_of_images):
        evt=run.event(times[t])
        img = xtcav_camera.image(evt)
        # skip if empty image
        if img is None: 
            continue

        xN = roiXN(evt)  #Size of the image in X                           
        x0 = roiX(evt)    #Position of the first pixel in x
        yN = roiYN(evt)  #Size of the image in Y 
        y0 = roiY(evt)    #Position of the first pixel in y
        
        if xN is None:       
            warnings.warn_explicit('No XTCAV ROI info',UserWarning,'XTCAV',0)
            continue
            
        x = x0+np.arange(0, xN) 
        y = y0+np.arange(0, yN) 

        ROI_XTCAV = ROIMetrics(xN, x0, yN, y0, x, y, valid=1) 

        return ROI_XTCAV, t

    return ROIMetrics(), end_of_images 


def GetShotToShotParameters(ebeam, gasdetector, evt_id):
    time = evt_id.time()
    sec  = time[0]
    nsec = time[1]
    unixtime = int((sec<<32)|nsec)
    fiducial = evt_id.fiducials()
 
    shot_to_shot = ShotToShotParameters(unixtime = unixtime, fiducial = fiducial)

    if ebeam:    
        ebeamcharge=ebeam.ebeamCharge()
        xtcavrfamp=ebeam.ebeamXTCAVAmpl()
        xtcavrfphase=ebeam.ebeamXTCAVPhase()
        dumpecharge=ebeam.ebeamDumpCharge()*Constants.E_CHARGE #In C  
        shot_to_shot = shot_to_shot._replace(ebeamcharge = ebeamcharge, 
            xtcavrfphase = xtcavrfphase, xtcavrfamp = xtcavrfamp, dumpecharge = dumpecharge)      
    else:    
        warnings.warn_explicit('No ebeamv info',UserWarning,'XTCAV',0)
        shot_to_shot = shot_to_shot._replace(valid = 0)
        
    if gasdetector:
        energydetector=(gasdetector.f_11_ENRC()+gasdetector.f_12_ENRC())/2    
    else:   #Some hardcoded values
        energydetector = Constants.ENERGY_DETECTOR
        warnings.warn_explicit('No gas detector info',UserWarning,'XTCAV',0)
        shot_to_shot = shot_to_shot._replace(valid = 0) 

    shot_to_shot = shot_to_shot._replace(xrayenergy = 1e-3*energydetector)

    return shot_to_shot


def DivideImageTasks(first_image, last_image, rank, size):
    num_shots = last_image - first_image
    tiling = np.arange(rank*4, rank*4+4,1) #  returns [0, 1, 2, 3] if e.g. rank == 0 and size == 4:
    comb1 = np.tile(tiling, np.ceil(num_shots/(4.*size)).astype(int))  # returns [0, 1, 2, 3, 0, 1, 2, 3, ...]        
    comb2 = np.repeat(np.arange(0, np.ceil(num_shots/(4.*size)), 1), 4) # returns [0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, ...]
            #  list of shot numbers assigned to this core
    main = comb2*4*size + comb1  + first_image # returns [  0.   1.   2.   3.  16.  17.  18.  19.  32.  33. ... ]
    main = np.delete(main, np.where(main>=last_image) )  # remove element if greater or equal to maximum number of shots in run
    return main.astype(int)
    
