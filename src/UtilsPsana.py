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
       
        globalCalibration = GlobalCalibration(
            umperpix=umperpix(evt), 
            strstrength=strstrength(evt), 
            rfampcalib=rfampcalib(evt), 
            rfphasecalib=rfphasecalib(evt), 
            dumpe=dumpe(evt), 
            dumpdisp=dumpdisp(evt) 
        ) 
        if globalCalibration.valid:
            return globalCalibration, t
                
    return globalCalibration, -1
          

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
       
        ROI_XTCAV = ROIMetrics(
            roiXN(evt), 
            roiX(evt), 
            roiYN(evt), 
            roiY(evt)
        ) 

        if ROI_XTCAV.valid: 
            return ROI_XTCAV, t

    return ROI_XTCAV, -1
    
def ShotToShotParameters(ebeam,gasdetector):
    """
    Obtain shot to shot parameters
    Arguments:
      ebeam: psana.Bld.BldDataEBeamV5
      gasdetector: psana.Bld.BldDataFEEGasDetEnergy
    Output:
      ROI: struct with the ROI
      ok: if all the data was retrieved correctly
    """
    ok=1
    echarge=1.60217657e-19;

    if ebeam:    
        ebeamcharge=ebeam.ebeamCharge()
        xtcavrfamp=ebeam.ebeamXTCAVAmpl()
        xtcavrfphase=ebeam.ebeamXTCAVPhase()
        dumpecharge=ebeam.ebeamDumpCharge()*echarge #In C        
    else:    #Some hardcoded values
        warnings.warn_explicit('No ebeamv info',UserWarning,'XTCAV',0)
        ok=0
        ebeamcharge=5
        xtcavrfamp=20
        xtcavrfphase=90
        energydetector=0.2;
        dumpecharge=175e-12 #In C
        
    if gasdetector:
        energydetector=(gasdetector.f_11_ENRC()+gasdetector.f_12_ENRC())/2    
    else:   #Some hardcoded values
        warnings.warn_explicit('No gas detector info',UserWarning,'XTCAV',0)
        ok=0
        energydetector=0.2;
        
    energy=1e-3*energydetector #In J
            
    shotToShot={
        'ebeamcharge':ebeamcharge,  #ebeamcharge
        'dumpecharge':dumpecharge,  #dumpecharge in C
        'xtcavrfamp': xtcavrfamp,   #RF amplitude
        'xtcavrfphase':xtcavrfphase, #RF phase
        'xrayenergy':energy         #Xrays energy in J
        }        
       
    return shotToShot,ok 

def DivideImageTasks(num_shots, rank, size):
    tiling = np.arange(rank*4, rank*4+4,1) #  returns [0, 1, 2, 3] if e.g. rank == 0 and size == 4:
    comb1 = np.tile(tiling, np.ceil(num_shots/(4.*size)).astype(np.int))  # returns [0, 1, 2, 3, 0, 1, 2, 3, ...]        
    comb2 = np.repeat(np.arange(0, np.ceil(num_shots/(4.*size)), 1), 4) # returns [0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, ...]
            #  list of shot numbers assigned to this core
    main = comb2*4*size + comb1  # returns [  0.   1.   2.   3.  16.  17.  18.  19.  32.  33. ... ]
    main = np.delete(main, np.where(main>=num_shots))  # remove element if greater or equal to maximum number of shots in run
    return main
    
