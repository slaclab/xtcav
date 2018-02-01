#(c) Coded by Alvaro Sanchez-Gonzalez 2014
import os
import time
import psana
import numpy as np
import glob
import sys
import getopt
import warnings
import Utils as xtu
import UtilsPsana as xtup
from DarkBackground import *
from CalibrationPaths import *

"""
    Function that generates a dark background image for XTCAV reconstruction purposes
    Arguments:
        experiment (str): String with the experiment reference to use. E.g. 'amoc8114'
        run (str): String with a run number. E.g. '123' 
        maxshots (int): Maximum number of images to use for the reference.
        calibrationpath (str): Custom calibration directory in case the default is not intended to be used.
        validityrange (tuple): If not set, the validity range for the reference will go from the 
        first run number used to generate the reference and the last run.
"""
def generateDarkBackground(
    experiment='amoc8114', 
    maxshots=401, 
    run_number='86', 
    validityrange=None, 
    calibrationpath='', 
    savetofile=True):
    
    warnings.filterwarnings('always',module='Utils',category=UserWarning)
    warnings.filterwarnings('ignore',module='Utils',category=RuntimeWarning, message="invalid value encountered in divide")
    
    """
    After setting all the parameters, this method has to be called to generate the dark reference and 
    save it in the proper location. 
    """
    print 'dark background reference'
    print '\t Experiment: %s' % experiment
    print '\t Run: %s' % run_number
    print '\t Valid shots to process: %d' % maxshots
    
    #Loading the dataset from the "dark" run, this way of working should be compatible with both xtc and hdf5 files
    dataSource=psana.DataSource("exp=%s:run=%s:idx" % (experiment, run_number))
    
    #Camera and type for the xtcav images
    xtcav_camera = psana.Source('DetInfo(XrayTransportDiagnostic.0:Opal1000.0)')
    xtcav_type=psana.Camera.FrameV1
    #xtcav_camera = psana.Detector('XrayTransportDiagnostic.0:Opal1000.0')
    
    #Stores for environment variables    
    configStore=dataSource.env().configStore();
    epicsStore=dataSource.env().epicsStore();

    db=DarkBackground()

    # ROI_XTCAV,ok=xtup.GetXTCAVImageROI(epicsStore)             
    # if not ok: #If the information is not good, we try next event
    #     return db
    # accumulator_xtcav=np.zeros(( ROI_XTCAV['yN'],ROI_XTCAV['xN']), dtype=np.float64)

    n=0  #Counter for the total number of xtcav images processed
     
    run = dataSource.runs().next()        
        #for e, evt in enumerate(dataSource.events()):
    times = run.times()
    for t in range(len(times)-1,-1,-1): #Starting from the back, to avoid waits in the cases where there are not xtcav images for the first shots
        evt=run.event(times[t])
    
        #ignore shots without xtcav, because we can get
        #incorrect EPICS information (e.g. ROI).  this is
        #a workaround for the fact that xtcav only records
        #epics on shots where it has camera data, as well
        #as an incorrect design in psana where epics information
        #is not stored per-shot (it is in a more global object
        #called "Env")
        #frame = xtcav_camera.image(evt)
        frame = evt.get(xtcav_type, xtcav_camera) 
        if frame is None: 
            continue

        if not 'ROI_XTCAV' in locals():   #After the first event the epics store should contain the ROI of the xtcav images, that let us get the x and y vectors
            ROI_XTCAV,ok=xtup.GetXTCAVImageROI(epicsStore)             
            if not ok: #If the information is not good, we try next event
                del ROI_XTCAV
                continue
            accumulator_xtcav=np.zeros(( ROI_XTCAV['yN'],ROI_XTCAV['xN']), dtype=np.float64)
                    
        #For each shot that contains an xtcav frame we retrieve it and add it to the accumulators
        img=frame.data16().astype(np.float64)
        
        n += 1
        accumulator_xtcav = accumulator_xtcav+img 
            
        if n % 5 == 0:
            sys.stdout.write('\r%.1f %% done, %d / %d' % ( float(n) / maxshots*100, n, maxshots ))
            sys.stdout.flush()   
        if n>=maxshots:                    #After a certain number of shots we stop (Ideally this would be an argument, rather than a hardcoded value)
            sys.stdout.write('\n')
            break                          
    #At the end of the program the total accumulator is saved     
    db.n=n
    db.image=accumulator_xtcav/n
    db.ROI=ROI_XTCAV
    db.run=run_number
    
    if not validityrange:
        validityrange=[run_number, 'end']
        
    cp=CalibrationPaths(dataSource.env(),calibrationpath)
    file=cp.newCalFileName('pedestals',validityrange[0],validityrange[1])
    
    if savetofile:
        db.Save(file)
    return db
        
    
    