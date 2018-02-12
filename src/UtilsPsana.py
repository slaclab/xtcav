import numpy as np
import psana
import warnings
from Metrics import *

def GetCameraSaturationValue(epicsStore, run, xtcav_camera, start=None):
    analysis_version = psana.Detector('XTCAV_Analysis_Version')
    times = run.times()

    if not start in range(-1, len(times)):
        start = len(times) - 1

    for t in range(start,-1,-1):
        evt=run.event(times[t])
        img = xtcav_camera.image(evt)
        # skip if empty image
        if img is None: 
            continue

        if analysis_version(evt) is not None:
            return (1<<12)-1
        else:
            return (1<<14)-1

def GetGlobalXTCAVCalibration(epicsStore, run, xtcav_camera, start=None):
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

    if not start in range(-1, len(times)):
        start = len(times) - 1

    for t in range(start,-1,-1):
        evt=run.event(times[t])
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
        valid = 1
        for k,v in global_calibration._asdict().iteritems():
            if not v:
                warnings.warn_explicit('No XTCAV Calibration for epics variable' + k,UserWarning,'XTCAV',0)
                valid = 0
        global_calibration._replace(valid=valid)

        if global_calibration.valid:
            return global_calibration, t
                
    return global_calibration, -1
          

def GetXTCAVImageROI(epicsStore, run, xtcav_camera, start=None):

    roiXN=psana.Detector('XTCAV_ROI_sizeX')
    roiX=psana.Detector('XTCAV_ROI_startX')
    roiYN=psana.Detector('XTCAV_ROI_sizeY')
    roiY=psana.Detector('XTCAV_ROI_startY')
    times = run.times()

    if not start in range(-1, len(times)):
        start = len(times) - 1

    for t in range(start,-1,-1):
        evt=run.event(times[t])
        img = xtcav_camera.image(evt)
        # skip if empty image
        if img is None: 
            continue

        xN = roiXN(evt)  #Size of the image in X                           
        x0 = roiX(evt)    #Position of the first pixel in x
        yN = roiYN(evt)  #Size of the image in Y 
        y0 = roiY(evt)    #Position of the first pixel in y
        
        valid = 1
        if roiX is None:       
            warnings.warn_explicit('No XTCAV ROI info',UserWarning,'XTCAV',0)
            valid = 0
            xN = 1024   #Size of the image in X                           
            x0 = 0         #Position of the first pixel in x
            yN = 1024   #Size of the image in Y 
            y0 = 0         #Position of the first pixel in y

        x = x0+np.arange(0, xN) 
        y = y0+np.arange(0, yN) 

        ROI_XTCAV = ROIMetrics(xN, x0, yN, y0, x, y, valid) 

        if valid: 
            return ROI_XTCAV, t

    return ROI_XTCAV, -1

def GetShotToShotParameters(ebeam, gasdetector, evt_id):
    ### move to constants file
    echarge=1.60217657e-19

    #Some default values
    ### move to constants file
    ebeamcharge=5
    xtcavrfamp=20
    xtcavrfphase=90
    energydetector=0.2
    dumpecharge=175e-12 #In C
    energydetector=0.2
    valid = 1

    if ebeam:    
        ebeamcharge=ebeam.ebeamCharge()
        xtcavrfamp=ebeam.ebeamXTCAVAmpl()
        xtcavrfphase=ebeam.ebeamXTCAVPhase()
        dumpecharge=ebeam.ebeamDumpCharge()*echarge #In C        
    else:    
        warnings.warn_explicit('No ebeamv info',UserWarning,'XTCAV',0)
        valid=0
        
    if gasdetector:
        energydetector=(gasdetector.f_11_ENRC()+gasdetector.f_12_ENRC())/2    
    else:   #Some hardcoded values
        warnings.warn_explicit('No gas detector info',UserWarning,'XTCAV',0)
        valid=0     

    xrayenergy=1e-3*energydetector #In J

    time = evt_id.time()
    sec  = time[0]
    nsec = time[1]
    unixtime = int((sec<<32)|nsec)
    fiducial = evt_id.fiducials()

    return ShotToShotParameters(
        ebeamcharge, dumpecharge, xtcavrfamp, 
        xtcavrfphase, xrayenergy,
        unixtime, fiducial, valid)

def DivideImageTasks(num_shots, rank, size):
    tiling = np.arange(rank*4, rank*4+4,1) #  returns [0, 1, 2, 3] if e.g. rank == 0 and size == 4:
    comb1 = np.tile(tiling, np.ceil(num_shots/(4.*size)).astype(np.int))  # returns [0, 1, 2, 3, 0, 1, 2, 3, ...]        
    comb2 = np.repeat(np.arange(0, np.ceil(num_shots/(4.*size)), 1), 4) # returns [0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, ...]
            #  list of shot numbers assigned to this core
    main = comb2*4*size + comb1  # returns [  0.   1.   2.   3.  16.  17.  18.  19.  32.  33. ... ]
    main = np.delete(main, np.where(main>=num_shots))  # remove element if greater or equal to maximum number of shots in run
    return main
    
