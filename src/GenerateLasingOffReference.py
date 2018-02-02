#(c) Coded by Alvaro Sanchez-Gonzalez 2014
# revised 31/07/15 by andr0s & polo5 to include parallel processing
import os
import time
import psana
import numpy as np
import glob
import pdb
import IPython
import sys
import getopt
import warnings
import Utils as xtu
import UtilsPsana as xtup
from DarkBackground import *
from LasingOffReference import *
from CalibrationPaths import *

# PP imports
from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
#print 'Core %s ... ready' % (rank + 1) # useful for debugging purposes
#sys.stdout.flush()

"""
    Function that generates a set of lasing off references for XTCAV reconstruction purposes
    Attributes:
        experiment (str): String with the experiment reference to use. E.g. 'amoc8114'
        runs (str): String with a run number, or a run interval. E.g. '123'  '134-156' 145,136'
        maxshots (int): Maximum number of images to use for the references.
        calibrationpath (str): Custom calibration directory in case the default is not intended to be used.
        nb (int): Number of bunches.
        medianfilter (int): Number of neighbours for median filter.
        snrfilter (float): Number of sigmas for the noise threshold.
        groupsize (int): Number of profiles to average together for each reference.
        roiwaistthres (float): ratio with respect to the maximum to decide on the waist of the XTCAV trace.
        roiexpand (float): number of waists that the region of interest around will span around the center of the trace.
        islandsplitmethod (str): island splitting algorithm. Set to 'scipylabel' or 'contourLabel'  The defaults parameter is 'scipylabel'.
"""

def generateLasingOffReference(
        experiment='amoc8114',             #Experiment label
        maxshots=401,                      #Maximum number of valid shots to process
        run_number='86',                         #Runs
        validityrange=None,
        darkreferencepath=None,             #Dark reference information
        nb=1,                              #Number of bunches
        groupsize=5,                       #Number of profiles to average together
        medianfilter=3,                    #Number of neighbours for median filter
        snrfilter=10,                      #Number of sigmas for the noise threshold
        roiwaistthres=0.2,                 #Parameter for the roi location
        roiexpand=2.5,                     #Parameter for the roi location
        islandsplitmethod = 'scipyLabel',  #Method for island splitting
        islandsplitpar1 = 3.0,                      #Ratio between number of pixels between largest and second largest groups when calling scipy.label
        islandsplitpar2 = 5.0,                      #Ratio between number of pixels between second/third largest groups when calling scipy.label
        calpath='',
        savetofile=True

    ):
    
    #Handle warnings
    warnings.filterwarnings('always',module='Utils',category=UserWarning)
    warnings.filterwarnings('ignore',module='Utils',category=RuntimeWarning, message="invalid value encountered in divide")

    #Some default values for the options
    
    ROI_XTCAV = None
    global_calibration = None
    dark_background = None
    
    print 'Lasing off reference'
    print '\t Experiment: %s' % experiment
    print '\t Runs: %s' % run_number
    print '\t Number of bunches: %d' % nb
    print '\t Valid shots to process: %d' % maxshots
    print '\t Dark reference run: %s' % darkreferencepath
    
    #Loading the data, this way of working should be compatible with both xtc and hdf5 files
    dataSource=psana.DataSource("exp=%s:run=%s:idx" % (experiment, run_number))

    #Camera and type for the xtcav images
    xtcav_camera = psana.Source('DetInfo(XrayTransportDiagnostic.0:Opal1000.0)')
    #xtcav_camera = psana.Detector('XrayTransportDiagnostic.0:Opal1000.0')
    #xtcav_camera.set_do_offset(do_offset=False)
    xtcav_type=psana.Camera.FrameV1

    #Ebeam type: it should actually be the version 5 which is the one that contains xtcav stuff
    ebeam_data = psana.Detector('EBeam')

    #Gas detectors for the pulse energies
    gasdetector_data = psana.Detector('FEEGasDetEnergy')

    #Stores for environment variables   
    epicsStore = dataSource.env().epicsStore();

    n=0 #Counter for the total number of xtcav images processed

    #Empty lists for the statistics obtained from each image, the shot to shot properties, and the ROI of each image (although this ROI is initially the same for each shot, it becomes different when the image is cropped around the trace)
    listImageStats=[];
    listShotToShot=[];
    listROI=[];
    listPU=[]

        
    run=dataSource.runs().next(); #This line and the previous line are a temporal hack to go only through the first run, that avoids an unexpected block when calling next at the iterator, when there are not remaining runs.
    num_processed = 0 #Counter for the total number of xtcav images processed within the run        
    times = run.times()

    #  Parallel Processing implementation by andr0s and polo5
    #  The run will be segmented into chunks of 4 shots, with each core alternatingly assigned to each.
    #  e.g. Core 1 | Core 2 | Core 3 | Core 1 | Core 2 | Core 3 | ....
    num_shots = len(times) #  The number of shots in this run
    image_numbers = xtup.DivideImageTasks(num_shots, rank, size)
    current_shot = 0

    for t in image_numbers[::-1]: #  Starting from the back, to avoid waits in the cases where there are not xtcav images for the first shots
        evt=run.event(times[int(t)])

        #ignore shots without xtcav, because we can get
        #incorrect EPICS information (e.g. ROI).  this is
        #a workaround for the fact that xtcav only records
        #epics on shots where it has camera data, as well
        #as an incorrect design in psana where epics information
        #is not stored per-shot (it is in a more global object
        #called "Env")
        frame = evt.get(xtcav_type, xtcav_camera) 
        if frame is None: 
            continue

        ebeam = ebeam_data.get(evt)
        gasdetector = gasdetector_data.get(evt)

        if not ROI_XTCAV: 
            roi, ok = xtup.GetXTCAVImageROI(epicsStore) 
            if not ok: #If the information is not good, we try next event
                continue
            ROI_XTCAV = roi

        if not global_calibration:
            gl_cal, ok = xtup.GetGlobalXTCAVCalibration(epicsStore)
            if not ok: #If the information is not good, we try next event
                continue
            global_calibration = gl_cal
            saturationValue = xtup.GetCameraSaturationValue(epicsStore)

        #If we have not loaded the dark background information yet, we do
        if dark_background is None:
            if not darkreferencepath:
                cp = CalibrationPaths(dataSource.env(), calpath)
                darkreferencepath = cp.findCalFileName('pedestals',evt.run())
                
            if not darkreferencepath:
                print ('Dark reference for run %d not found, image will not be background substracted' % evt.run())
                loadeddarkreference=False 
                dark_background = False
            else:
                dark_background = DarkBackground.Load(darkreferencepath)       
                      
        if frame: #For each shot that contains an xtcav frame we retrieve it        
            img=frame.data16().astype(np.float64)

            if np.max(img)>=saturationValue : #Detection if the image is saturated, we skip if it is
                warnings.warn_explicit('Saturated Image',UserWarning,'XTCAV',0)
                continue

            shotToShot,ok = xtup.ShotToShotParameters(ebeam, gasdetector) #Obtain the shot to shot parameters necessary for the retrieval of the x and y axis in time and energy units
            if not ok: #If the information is not good, we skip the event
                continue

            id = evt.get(psana.EventId)
            time = id.time()
            sec  = time[0]
            nsec = time[1]
            shotToShot['unixtime'] = int((sec<<32)|nsec)
            shotToShot['fiducial'] = id.fiducials()
            
            #Subtract the dark background, taking into account properly possible different ROIs, if it is available
            if dark_background:        
                img, ROI=xtu.SubtractBackground(img, ROI_XTCAV, dark_background.image, dark_background.ROI) 
            else:
                ROI = ROI_XTCAV
            img, ok=xtu.DenoiseImage(img, medianfilter, snrfilter)                    #Remove noise from the image and normalize it
            if not ok:                                        #If there is nothing in the image we skip the event  
                continue

            img, ROI=xtu.FindROI(img, ROI, roiwaistthres, roiexpand)                  #Crop the image, the ROI struct is changed. It also add an extra dimension to the image so the array can store multiple images corresponding to different bunches
            if ROI['xN']<3 or ROI['yN']<3:
                print 'ROI too small',ROI['xN'],ROI['yN']
                continue

            img = xtu.SplitImage(img, nb, islandsplitmethod, islandsplitpar1, islandsplitpar2)#new

            if nb!=img.shape[0]:
                continue
            imageStats=xtu.ProcessXTCAVImage(img,ROI)          #Obtain the different properties and profiles from the trace               

            PU, ok = xtu.CalculatePhysicalUnits(ROI,[imageStats[0]['xCOM'],imageStats[0]['yCOM']],shotToShot,global_calibration)   
            if not ok:
                continue

            #If the step in time is negative, we mirror the x axis to make it ascending and consequently mirror the profiles
            if PU['xfsPerPix']<0:
                PU['xfs']=PU['xfs'][::-1]
                NB=len(imageStats)
                for j in range(NB):
                    imageStats[j]['xProfile']=imageStats[j]['xProfile'][::-1]
                    imageStats[j]['yCOMslice']=imageStats[j]['yCOMslice'][::-1]
                    imageStats[j]['yRMSslice']=imageStats[j]['yRMSslice'][::-1]                                               
                                                                                                                                                                                    
            listImageStats.append(imageStats)
            listShotToShot.append(shotToShot)
            listROI.append(ROI)
            listPU.append(PU)
            
            n += 1
            num_processed += 1
            # print core numb and percentage

            if current_shot % 5 == 0:
                if size==1:extrainfo='\r'
                else:extrainfo='\nCore %d: '%(rank + 1)
                sys.stdout.write('%s%.1f %% done, %d / %d' % (extrainfo, float(current_shot) / np.ceil(maxshots/float(size)) *100, current_shot, np.ceil(maxshots/float(size))))
                sys.stdout.flush()
            current_shot += 1
            if current_shot >= np.ceil(maxshots/float(size)):
                sys.stdout.write('\n')
                break

    #  here gather all shots in one core, add all lists
    exp = {'listImageStats': listImageStats, 'listShotToShot': listShotToShot, 'listROI': listROI, 'listPU': listPU}
    processedlist = comm.gather(exp, root=0)
    
    if rank != 0:
        return
    
    listImageStats = []
    listShotToShot = []
    listROI = []
    listPU = []
    
    for i in range(size):
        p = processedlist[i]
        listImageStats += p['listImageStats']
        listShotToShot += p['listShotToShot']
        listROI += p['listROI']
        listPU += p['listPU']
        
    #Since there are 12 cores it is possible that there are more references than needed. In that case we discard some
    n=len(listImageStats)
    if n>maxshots:
        n=maxshots
        listImageStats=listImageStats[0:n]
        listShotToShot=listShotToShot[0:n]
        listROI=listROI[0:n]
        listPU=listPU[0:n]
        
    #At the end, all the reference profiles are converted to Physical units, grouped and averaged together
    averagedProfiles = xtu.AverageXTCAVProfilesGroups(listROI,listImageStats,listShotToShot,listPU,groupsize);     

    lor=LasingOffReference()
    lor.averagedProfiles=averagedProfiles
    lor.runs=runs
    lor.n=n
    
    # n should be consistent with len(final list)
    
    parameters= {
        'version' : 0,
        'darkreferencepath':darkreferencepath,
        'nb':nb,
        'groupsize':groupsize,
        'medianfilter':medianfilter,
        'snrfilter':snrfilter,
        'roiwaistthres':roiwaistthres,
        'roiexpand':roiexpand,
        'islandsplitmethod':islandsplitmethod,
        'islandsplitpar1':islandsplitpar1,
        'islandsplitpar2':slandsplitpar2,
    }
    
    
    lor.parameters=parameters
    if not validityrange:
        validityrange=[runs[0], 'end']
        
    cp=CalibrationPaths(dataSource.env(),calpath)
    file=cp.newCalFileName('lasingoffreference',validityrange[0],validityrange[1])
       
    if savetofile:
        lor.Save(file)
        
   
    #Access to the different properties: here we can change flags when a parameter is changed, or check the validity of the property

