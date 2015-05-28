#(c) Coded by Alvaro Sanchez-Gonzalez 2014


import os
import time
import psana
import numpy as np
import glob
import sys
import getopt
import warnings
import Utils as xtu   #HACK
from DarkBackground import *
from CalibrationPaths import *


class GenerateDarkBackground(object):
    """
    Class that generates a dark background image for XTCAV reconstruction purposes
    Attributes:
        experiment (str): String with the experiment reference to use. E.g. 'amoc8114'
        runs (str): String with a run number, or a run interval. E.g. '123'  '134-156' 145,136'
        maxshots (int): Maximum number of images to use for the reference.
        calibrationpath (str): Custom calibration directory in case the default is not intended to be used.
    """

    def __init__(self):
    
        #Handle warnings
        warnings.filterwarnings('always',module='Utils',category=UserWarning)
        warnings.filterwarnings('ignore',module='Utils',category=RuntimeWarning, message="invalid value encountered in divide")
        
        self._experiment='amoc8114'
        self._maxshots=100
        self._runs='85'
        self._validityrange=[]
        self._calpath=''
                       
    def Generate(self):
        """
        After setting all the parameters, this method has to be called to generate the dark reference and save it in the proper location. It not set, the validity range for the reference will go from the first run number used to generate the reference and the last run.
        """
        print 'dark background reference'
        print '\t Experiment: %s' % self._experiment
        print '\t Runs: %s' % self._runs
        print '\t Valid shots to process: %d' % self._maxshots
        
        #Loading the dataset from the "dark" run, this way of working should be compatible with both xtc and hdf5 files
        dataSource=psana.DataSource("exp=%s:run=%s:idx" % (self._experiment,self._runs))
        
        #Camera and type for the xtcav images
        xtcav_camera = psana.Source('DetInfo(XrayTransportDiagnostic.0:Opal1000.0)')
        xtcav_type=psana.Camera.FrameV1
        
        #Stores for environment variables    
        configStore=dataSource.env().configStore();
        epicsStore=dataSource.env().epicsStore();

        n=0  #Counter for the total number of xtcav images processed
         
        runs=numpy.array([],dtype=int) #Array that contains the run processed run numbers

        #for r,run in enumerate(dataSource.runs()): 
        for r in [0]:
            run=dataSource.runs().next(); #This line and the previous line are a temporal hack to go only through the first run, that avoids an unexpected block when calling next at the iterator, when there are not remaining runs.
            runs = numpy.append(runs,run.run());
            n_r=0        #Counter for the total number of xtcav images processed within the run
            
            #for e, evt in enumerate(dataSource.events()):
            times = run.times()
            for t in range(len(times)-1,-1,-1): #Starting from the back, to avoid waits in the cases where there are not xtcav images for the first shots
                evt=run.event(times[t])
            
                if not 'ROI_XTCAV' in locals():   #After the first event the epics store should contain the ROI of the xtcav images, that let us get the x and y vectors
                    ROI_XTCAV,ok=xtu.GetXTCAVImageROI(epicsStore)             
                    if not ok: #If the information is not good, we try next event
                        del ROI_XTCAV
                        continue
                    accumulator_xtcav=np.zeros(( ROI_XTCAV['yN'],ROI_XTCAV['xN']), dtype=np.float64)
                            
                frame = evt.get(xtcav_type, xtcav_camera) 
                if frame:                       #For each shot that contains an xtcav frame we retrieve it and add it to the accumulators
                    img=frame.data16().astype(np.float64)
                    
                    n=n+1
                    n_r=n_r+1
                    print "%d/%d" % (n_r,self._maxshots) #Debugging purposes, will be removed
                    
                    accumulator_xtcav=accumulator_xtcav+img 
                                           
                if n_r>=self._maxshots:                    #After a certain number of shots we stop (Ideally this would be an argument, rather than a hardcoded value)
                    break      

        
        #At the end of the program the total accumulator is saved     
        db=DarkBackground()
        db.n=n
        db.image=accumulator_xtcav/n
        db.ROI=ROI_XTCAV
        db.runs=runs
        
        if not self._validityrange:
            validityrange=[runs[0], 'end']
        else:
            validityrange=self._validityrange
            
        cp=CalibrationPaths(dataSource.env(),self._calpath)
        file=cp.newCalFileName('pedestals',validityrange[0],validityrange[1])
        
        db.Save(file)
        
    def SetValidityRange(self,runBegin,runEnd='end'):
        """Sets the validity range for the generated reference.

        Args:
            runBegin (int): First run in the range.
            runEnd (int or str): Last run in the range (use 'end' to leave the range open)..

        """
        self._validityrange=[runBegin, runEnd];
    
    #Access to the different properties: here we can change flags when a parameter is changed, or check the validity of the property
    @property
    def maxshots(self):
        return self._maxshots
    @maxshots.setter
    def maxshots(self, maxshots):
        self._maxshots = maxshots     
    @property        
    def experiment(self):
        return self._experiment
    @experiment.setter
    def experiment(self, experiment):
        self._experiment = experiment
    @property
    def runs(self):
        return self._runs
    @runs.setter
    def runs(self, runs):
        self._runs = runs       
    @property
    def calibrationpath(self):
        return self._calpath
    @calibrationpath.setter
    def calibrationpath(self, calpath):
        self._calpath = calpath

