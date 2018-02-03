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
    print "here"
    
    #Camera and type for the xtcav images
    xtcav_camera = psana.Detector('XrayTransportDiagnostic.0:Opal1000.0')
    
    #Stores for environment variables    
    configStore=dataSource.env().configStore();
    epicsStore=dataSource.env().epicsStore();

    n=0  #Counter for the total number of xtcav images processed 
    run = dataSource.runs().next()        
    

    ROI_XTCAV, last_image = xtup.GetXTCAVImageROI(epicsStore, run, xtcav_camera)
    accumulator_xtcav=np.zeros((ROI_XTCAV.YN, ROI_XTCAV.XN), dtype=np.float64)

    times = run.times()
    for t in range(last_image,-1,-1): #Starting from the last valid image, to avoid waits in the cases where there are not xtcav images for the first shots
        evt=run.event(times[t])
    
        #ignore shots without xtcav, because we can get incorrect EPICS information (e.g. ROI).  this is
        #a workaround for the fact that xtcav only records epics on shots where it has camera data, as well
        #as an incorrect design in psana where epics information is not stored per-shot (it is in a more global object
        #called "Env")
        img = xtcav_camera.image(evt)
        # skip if empty image
        if img is None: 
            continue
      
        accumulator_xtcav += img 
        n += 1
            
        if n % 5 == 0:
            sys.stdout.write('\r%.1f %% done, %d / %d' % ( float(n) / maxshots*100, n, maxshots ))
            sys.stdout.flush()   
        if n>=maxshots:                    #After a certain number of shots we stop (Ideally this would be an argument, rather than a hardcoded value)
            sys.stdout.write('\n')
            break                          
    #At the end of the program the total accumulator is saved  
    db=DarkBackground()   
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
        
    
    