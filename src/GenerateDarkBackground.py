#(c) Coded by Alvaro Sanchez-Gonzalez 2014


import os
import time
import psana
import numpy as np
import glob
import sys
import getopt
import xtcav.Utils as xtu
from xtcav.DarkBackground import *
from xtcav.CalibrationPaths import *


class GenerateDarkBackground(object):
    def __init__(self):
        self._experiment='amoc8114'
        self._maxshots=100
        self._runs='85'
        self._validityrange=[]
                       
    def Generate(self):
        print 'dark background reference'
        print '\t Experiment: %s' % self._experiment
        print '\t Runs: %s' % self._runs
        print '\t Valid shots to process: %d' % self._maxshots
        
        #Loading the dataset from the "dark" run, this way of working should be compatible with both xtc and hdf5 files
        dataSource=psana.DataSource("exp=%s:run=%s" % (self._experiment,self._runs))
        
        #Camera and type for the xtcav images
        xtcav_camera = psana.Source('DetInfo(XrayTransportDiagnostic.0:Opal1000.0)')
        xtcav_type=psana.Camera.FrameV1
        
        #Stores for environment variables    
        configStore=dataSource.env().configStore();
        epicsStore=dataSource.env().epicsStore();

        #Ebeam type: it should actually be the version 5 which is the one that contains xtcav stuff
        ebeam_data=psana.Source('BldInfo(EBeam)')
        ebeam_type=psana.Bld.BldDataEBeamV4

        n=0  #Counter for the total number of xtcav images processed
         
        runs=numpy.array([],dtype=int) #Array that contains the run processed run numbers

        #for r,run in enumerate(dataSource.runs()): 
        for r in [0]:
            run=dataSource.runs().next(); #This line and the previous line are a temporal hack to go only through the first run, that avoids an unexpected block when calling next at the iterator, when there are not remaining runs.
            runs = numpy.append(runs,run.run());
            n_r=0        #Counter for the total number of xtcav images processed within the run
            
            for e, evt in enumerate(dataSource.events()):
            
                ebeam = evt.get(ebeam_type,ebeam_data)    

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
            
        cp=CalibrationPaths(dataSource)
        file=cp.newCalFileName('pedestals',validityrange[0],validityrange[1])
        
        db.Save(file)
        
    def SetValidityRange(self,runBegin,runEnd='end'):
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

