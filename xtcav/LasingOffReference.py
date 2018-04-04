
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
import SplittingUtils as su
import Constants
from CalibrationPaths import *
from DarkBackground import *
from FileInterface import Load as constLoad
from FileInterface import Save as constSave

# PP imports
from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
#print 'Core %s ... ready' % (rank + 1) # useful for debugging purposes
#sys.stdout.flush()
"""
    Class that generates a set of lasing off references for XTCAV reconstruction purposes
    Attributes:
        experiment (str): String with the experiment reference to use. E.g. 'amoc8114'
        runs (str): String with a run number, or a run interval. E.g. '123'  '134-156' 145,136'
        maxshots (int): Maximum number of images to use for the references.
        start_image (int): image in run to start from
        validityrange (tuple): If not set, the validity range for the reference will go from the 
        first run number used to generate the reference and the last run.
        calibrationpath (str): Custom calibration directory in case the default is not intended to be used.
        num_bunches (int): Number of bunches.
        medianfilter (int): Number of neighbours for median filter.
        snrfilter (float): Number of sigmas for the noise threshold.
        groupsize (int): Number of profiles to average together for each reference.
        roiwaistthres (float): ratio with respect to the maximum to decide on the waist of the XTCAV trace.
        roiexpand (float): number of waists that the region of interest around will span around the center of the trace.
        islandsplitmethod (str): island splitting algorithm. Set to 'scipylabel' or 'contourLabel'  The defaults parameter is 'scipylabel'.
"""

class LasingOffReference(object):

    def __init__(self,
            experiment='amoc8114',  #Experiment label
            maxshots=401,           #Maximum number of valid shots to process
            run_number='86',        #Run number
            start_image=0,          #Starting image in run
            validityrange=None,
            darkreferencepath=None, #Dark reference information
            num_bunches=1,                   #Number of bunches
            groupsize=0 ,           #Number of profiles to average together
            medianfilter=3,         #Number of neighbours for median filter in algorithm
            snrfilter=10,           #Number of sigmas for the noise threshold
            roiwaistthres=0.2,      #Parameter for the roi location
            roiexpand=2.5,          #Parameter for the roi location
            islandsplitmethod = Constants.DEFAULT_SPLIT_METHOD,      #Method for island splitting
            islandsplitpar1 = 3.0,  #Ratio between number of pixels between largest and second largest groups when calling scipy.label
            islandsplitpar2 = 5.,   #Ratio between number of pixels between second/third largest groups when calling scipy.label
            calpath='',
            savetofile=True):

        self.parameters = LasingOffParameters(experiment = experiment,
            maxshots = maxshots, run = run_number, start = start_image, validityrange = validityrange, 
            darkreferencepath = darkreferencepath, num_bunches = num_bunches, groupsize=groupsize, 
            medianfilter=medianfilter, snrfilter=snrfilter, roiwaistthres=roiwaistthres,
            roiexpand = roiexpand, islandsplitmethod=islandsplitmethod, islandsplitpar2 = islandsplitpar2,
            islandsplitpar1=islandsplitpar1, calpath=calpath, version=1)


        warnings.filterwarnings('always',module='Utils',category=UserWarning)
        warnings.filterwarnings('ignore',module='Utils',category=RuntimeWarning, message="invalid value encountered in divide")
        
        if rank == 0:
            print 'Lasing off reference'
            print '\t Experiment: %s' % self.parameters.experiment
            print '\t Runs: %s' % self.parameters.run
            print '\t Number of bunches: %d' % self.parameters.num_bunches
            print '\t Valid shots to process: %d' % self.parameters.maxshots
            print '\t Dark reference run: %s' % self.parameters.darkreferencepath
        
        #Loading the data, this way of working should be compatible with both xtc and hdf5 files
        dataSource = psana.DataSource("exp=%s:run=%s:idx" % (self.parameters.experiment, self.parameters.run))

        #Camera for the xtcav images
        xtcav_camera = psana.Detector(Constants.SRC)

        #Ebeam type
        ebeam_data = psana.Detector(Constants.EBEAM)

        #Gas detectors for the pulse energies
        gasdetector_data = psana.Detector(Constants.GAS_DETECTOR)

        #Empty list for the statistics obtained from each image, the shot to shot properties, and the ROI of each image (although this ROI is initially the same for each shot, it becomes different when the image is cropped around the trace)
        list_image_profiles= []
            
        run = dataSource.runs().next()
        env = dataSource.env()

        dark_background = self._getDarkBackground(env)

        #Calibration values needed to process images. first_event is the index of the first event with valid data
        roi_xtcav, global_calibration, saturation_value, first_event = self._getCalibrationValues(run, xtcav_camera)

        #  Parallel Processing implementation by andr0s and polo5
        #  The run will be segmented into chunks of 4 shots, with each core alternatingly assigned to each.
        #  e.g. Core 1 | Core 2 | Core 3 | Core 1 | Core 2 | Core 3 | ....
        times = run.times()
        image_numbers = xtup.DivideImageTasks(first_event, len(times), rank, size)

        num_processed = 0 #Counter for the total number of xtcav images processed within the run  
        for t in image_numbers: 
            evt = run.event(times[t])
            ebeam = ebeam_data.get(evt)
            gasdetector = gasdetector_data.get(evt)

            shot_to_shot = xtup.GetShotToShotParameters(ebeam, gasdetector, evt.get(psana.EventId)) #Obtain the shot to shot parameters necessary for the retrieval of the x and y axis in time and energy units
        
            if not shot_to_shot.valid: #If the information is not good, we skip the event
                continue 

            img = xtcav_camera.image(evt)
            image_profile = xtu.processImage(img, self.parameters, dark_background, global_calibration, 
                                                    saturation_value, roi_xtcav, shot_to_shot)

            if not image_profile:
                continue
            
            #Append only image profile, omit processed image                                                                                                                                                              
            list_image_profiles.append(image_profile[0])     
            num_processed += 1

            self._printProgressStatements(num_processed)

            if num_processed >= np.ceil(self.parameters.maxshots/float(size)):
                break

        # here gather all shots in one core, add all lists
        image_profiles = comm.gather(list_image_profiles, root=0)
        
        if rank != 0:
            return

        sys.stdout.write('\n')
        # Flatten gathered arrays
        image_profiles = [item for sublist in image_profiles for item in sublist]

        #Since there are 12 cores it is possible that there are more references than needed. In that case we discard some
        if len(image_profiles) > self.parameters.maxshots:
            image_profiles = image_profiles[0:self.parameters.maxshots]
         
        #At the end, all the reference profiles are converted to Physical units, grouped and averaged together
        averaged_profiles = xtu.AverageXTCAVProfilesGroups(image_profiles, self.parameters.groupsize);     

        self.averaged_profiles=averaged_profiles
        self.n=num_processed    
        
        if not self.parameters.validityrange:
            self.parameters = self.parameters._replace(validityrange=(self.parameters.run, 'end'))

        if savetofile:
            cp = CalibrationPaths(env, self.parameters.calpath)
            file = cp.newCalFileName('lasingoffreference', self.parameters.validityrange[0], self.parameters.validityrange[1])
            self.save(file)


    def _printProgressStatements(self, num_processed):
        # print core numb and percentage
        if num_processed % 5 == 0:
            extrainfo = '\r' if size == 1 else '\nCore %d: '%(rank + 1)
            sys.stdout.write('%s%.1f %% done, %d / %d' % (extrainfo, float(num_processed) / np.ceil(self.parameters.maxshots/float(size)) *100, num_processed, np.ceil(self.parameters.maxshots/float(size))))
            sys.stdout.flush()


    def _getDarkBackground(self, env):
        """
        Internal method. Loads dark background reference
        """
        if not self.parameters.darkreferencepath:
            cp = CalibrationPaths(env, self.parameters.calpath)
            darkreferencepath = cp.findCalFileName('pedestals', int(self.parameters.run))
            if not darkreferencepath:
                print ('Dark reference for run %s not found, image will not be background substracted' % self.parameters.run)
                return None

            self.parameters = self.parameters._replace(darkreferencepath = darkreferencepath)
        return DarkBackground.load(self.parameters.darkreferencepath)


    @staticmethod
    def _getCalibrationValues(run, xtcav_camera):
        roi_xtcav, global_calibration, saturation_value = None, None, None
        times = run.times()

        end_of_images = len(times)
        for t in range(end_of_images):
            evt = run.event(times[t])
            img = xtcav_camera.image(evt)
            # skip if empty image
            if img is None: 
                continue

            roi_xtcav = xtup.GetXTCAVImageROI(evt)
            global_calibration = xtup.GetGlobalXTCAVCalibration(evt)
            saturation_value = xtup.GetCameraSaturationValue(evt)

            if not roi_xtcav or not global_calibration or not saturation_value:
                continue

            return roi_xtcav, global_calibration, saturation_value, t

        return roi_xtcav, global_calibration, saturation_value, end_of_images


    def save(self, path):

        ###Move this to file interface folder...
        instance = copy.deepcopy(self)
        instance.parameters = dict(vars(self.parameters))
        instance.averaged_profiles = dict(vars(self.averaged_profiles))
        constSave(instance,path)

    @staticmethod
    def load(path):
        lor = constLoad(path)
        try:
            lor.parameters = LasingOffParameters(**lor.parameters)
            lor.averaged_profiles = xtu.AveragedProfiles(**lor.averaged_profiles)
        except TypeError:
            print "Could not load Lasing Off Reference with path "+ path+". Try recreating lasing off " +\
            "reference to ensure compatability between versions"
            sys.exit(1)
        return lor


LasingOffParameters = xtu.namedtuple('LasingOffParameters', 
    ['experiment', 
    'maxshots', 
    'run', 
    'start',
    'validityrange', 
    'darkreferencepath', 
    'num_bunches', 
    'groupsize', 
    'medianfilter', 
    'snrfilter', 
    'roiwaistthres', 
    'roiexpand', 
    'islandsplitmethod',
    'islandsplitpar1', 
    'islandsplitpar2', 
    'calpath', 
    'version'])

