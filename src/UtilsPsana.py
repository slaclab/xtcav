import numpy as np
import psana
import warnings
from ROIMetrics import *

def GetGlobalCalibValue(epicsStore,names,ok):
    """
    Iterate through list of names to try to get calibration value
    Arguments:
      names: list of string-names
    Output:
      value of epics variable and flag: 1 if found, 0 if not found.
      if not found a default value of 0 is returned for the epics variable.
    """
    for n in names:
        val = epicsStore.value(n)
        if val is not None: return val
    warnings.warn_explicit('No XTCAV Calibration for epics variable'+name[0],UserWarning,'XTCAV',0)
    ok[0]=0 # notify caller that no value was found for a variable
    return 0

def GetCameraSaturationValue(epicsStore):
    if epicsStore.value('XTCAV_Analysis_Version') is not None:
        return (1<<12)-1
    else:
        return (1<<14)-1

def GetGlobalXTCAVCalibration(epicsStore):
    """
    Obtain the global XTCAV calibration form the epicsStore
    Arguments:
      epicsStore
    Output:
      globalCalibration: struct with the parameters
      ok: if all the data was retrieved correctly
    """
    

    ok = [1]
    umperpix=GetGlobalCalibValue(epicsStore,['XTCAV_calib_umPerPx','OTRS:DMP1:695:RESOLUTION'],ok)
    strstrength=GetGlobalCalibValue(epicsStore,['XTCAV_strength_par_S','Streak_Strength','OTRS:DMP1:695:TCAL_X'],ok)
    rfampcalib=GetGlobalCalibValue(epicsStore,['XTCAV_Amp_Des_calib_MV','XTCAV_Cal_Amp','SIOC:SYS0:ML01:AO214'],ok)
    rfphasecalib=GetGlobalCalibValue(epicsStore,['XTCAV_Phas_Des_calib_deg','XTCAV_Cal_Phase','SIOC:SYS0:ML01:AO215'],ok)
    dumpe=GetGlobalCalibValue(epicsStore,['XTCAV_Beam_energy_dump_GeV','Dump_Energy','REFS:DMP1:400:EDES'],ok)
    dumpdisp=GetGlobalCalibValue(epicsStore,['XTCAV_calib_disp_posToEnergy','Dump_Disp','SIOC:SYS0:ML01:AO216'],ok)

    globalCalibration={
        'umperpix':umperpix, #Pixel size of the XTCAV camera
        'strstrength':strstrength,  #Strength parameter
        'rfampcalib':rfampcalib,    #Calibration of the RF amplitude
        'rfphasecalib':rfphasecalib,    #Calibration of the RF phase
        'dumpe':dumpe,                  #Beam energy: dump config
        'dumpdisp':dumpdisp             #Vertical position to energy: dispersion
        }
                
                
    return globalCalibration,ok[0]
          

def GetXTCAVImageROI(epicsStore, run, xtcav_camera):
    roiXN=psana.Detector('XTCAV_ROI_sizeX')# epicsStore.value('XTCAV_ROI_sizeX')
    roiX=psana.Detector('XTCAV_ROI_startX')#epicsStore.value('XTCAV_ROI_startX')
    roiYN=psana.Detector('XTCAV_ROI_sizeY')#epicsStore.value('XTCAV_ROI_sizeY')
    roiY=psana.Detector('XTCAV_ROI_startY')#epicsStore.value('XTCAV_ROI_startY')
    times = run.times()

    for t in range(len(times)-1,-1,-1):
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
            break

    return ROI_XTCAV, t
    
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
    
