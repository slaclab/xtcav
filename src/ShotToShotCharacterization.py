#(c) Coded by Alvaro Sanchez-Gonzalez 2014

#Script for the retrieval of the pulses shot to shot

import os
import time
import psana
import numpy as np
import glob
import pdb
import IPython
import sys
import getopt
import math
import xtcav.Utils as xtu
from xtcav.DarkBackground import *
from xtcav.LasingOffReference import *
from xtcav.CalibrationPaths import *

class ShotToShotCharacterization(object):

    def __init__(self):
        #Some default values for the options
        self._experiment='amoc8114'     #Experiment label
        self._darkreferencepath=[];         #Dark reference information
        self._lasingoffreferencepath=[];         #Dark reference information
        self._darkreference=[]
        self._lasingoffreference=[]
        self._nb=float('nan')                    #Number of bunches
        self._medianfilter=float('nan')          #Number of neighbours for median filter
        self._snrfilter=float('nan')             #Number of sigmas for the noise threshold
        self._roiwaistthres=float('nan')         #Parameter for the roi location
        self._roiexpand=float('nan')             #Parameter for the roi location
        self._currentevent=[]
        self._eventresultsstep1=[]
        self._eventresultsstep2=[]
        self._datasource=[]
        self._globalcalibration=[]
        self._roixtcav=[]
        
        #Different flags
        self._loadeddarkreference=False
        self._loadedlasingoffreference=False
        self._currenteventavailable=False
        self._currenteventprocessedstep1=False   #Step one is pure processing of the trace, without lasing off referencing
        self._currenteventprocessedstep2=False   #Step two is the processing of the profiles with respect to the reference
        self._datasourceinfo=False      

        #Camera and type for the xtcav images
        self.xtcav_camera = psana.Source('DetInfo(XrayTransportDiagnostic.0:Opal1000.0)')
        self.xtcav_type=psana.Camera.FrameV1
        self._frame=[]
        
        #Ebeam type: it should actually be the version 5 which is the one that contains xtcav stuff
        self.ebeam_data=psana.Source('BldInfo(EBeam)')
        self.ebeam_type=psana.Bld.BldDataEBeamV5
        self._ebeam=[]

        #Gas detectors for the pulse energies
        self.gasdetector_data=psana.Source('BldInfo(FEEGasDetEnergy)')
        self.gasdetector_type=psana.Bld.BldDataFEEGasDetEnergy
        self._gasdetector=[]
        
    def LoadDarkReference(self):
        if not self._darkreferencepath:
            cp=CalibrationPaths(self._datasource)       
            self._darkreferencepath=cp.findCalFileName('pedestals',self._currentevent.run())
        
        self._darkreference=DarkBackground.Load(self._darkreferencepath)
        self._loadeddarkreference=True
        
    def LoadLasingOffReference(self):
        if not self._lasingoffreferencepath:
            cp=CalibrationPaths(self._datasource)       
            self._lasingoffreferencepath=cp.findCalFileName('lasingoffreference',self._currentevent.run())
        self._lasingoffreference=LasingOffReference.Load(self._lasingoffreferencepath)
        self._loadedlasingoffreference=True
                                
        #Only use the parameters if they have not been manually set, except for the number of bunches. That one is mandatory.
        self._nb=self._lasingoffreference.parameters['nb']
        if math.isnan(self._medianfilter):
            self._medianfilter=self._lasingoffreference.parameters['medianfilter']
        if math.isnan(self._snrfilter):
            self._snrfilter=self._lasingoffreference.parameters['snrfilter']
        if math.isnan(self._roiwaistthres):
            self._roiwaistthres=self._lasingoffreference.parameters['roiwaistthres']
        if math.isnan(self._roiexpand):
            self._roiexpand=self._lasingoffreference.parameters['roiexpand']
        if not self._darkreferencepath:
            self._darkreferencepath=self._lasingoffreference.parameters['darkreferencepath']
                           
    def SetCurrentEvent(self,evt):
    
        ebeam = evt.get(self.ebeam_type,self.ebeam_data)        
        gasdetector=evt.get(self.gasdetector_type,self.gasdetector_data) 
        frame = evt.get(self.xtcav_type, self.xtcav_camera) 
        
        self._currenteventprocessedstep1=False
        self._currenteventprocessedstep2=False
        self._eventresultsstep1=[]
        self._eventresultsstep2=[]
        
        if (ebeam and gasdetector and frame):
            self._ebeam=ebeam
            self._gasdetector=gasdetector
            self._frame=frame
            self._currentevent=evt           
            self._currenteventavailable=True
            return True
        else:
            self._currenteventavailable=False
            return False
        
    def SetDataSource(self,datasource):
        self._datasource=datasource
        self._datasourceinfo=False  

    def GetFullResults(self):
        if not self._currenteventprocessedstep2:
            self.ProcessShotStep2()
            
        return self._eventresultsstep2    
        
    def ProcessShotStep1(self):
        
        if not self._currenteventavailable:
            return False
           
        #It is important that this is open first so the experiment name is set properly (important for loading references)   
        if not self._datasourceinfo:
            self._experiment=self._datasource.env().experiment()
            epicsstore=self._datasource.env().epicsStore();
            self._globalCalibration,ok1=xtu.GetGlobalXTCAVCalibration(epicsstore)
            self._roixtcav,ok2=xtu.GetXTCAVImageROI(epicsstore) 
            if ok1 and ok2: #If the information is not good, we try next event
                self._datasourceinfo=True
            else:
                return False
        
        #It is important that the lasing off reference is open first, because it may reset the lasing off reference that needs to be loaded        
        if not self._loadedlasingoffreference:
            self.LoadLasingOffReference()
        
        if not self._loadeddarkreference:
            self.LoadDarkReference()
                                                    
        img=self._frame.data16().astype(np.float64)           
        if np.max(img)>=16383 : #Detection if the image is saturated, we skip if it is
            print 'Saturated Image'
            
        shotToShot,ok = xtu.ShotToShotParameters(self._ebeam,self._gasdetector) #Obtain the shot to shot parameters necessary for the retrieval of the x and y axis in time and energy units
        if not ok: #If the information is not good, we skip the event
            return False
                            
        img,ROI=xtu.SubtractBackground(img,self._roixtcav,self._darkreference.image,self._darkreference.ROI)  #Subtract the dark background, taking into account properly possible different ROIs
        img,ok=xtu.DenoiseImage(img,self._medianfilter,self._snrfilter)                    #Remove noise from the image and normalize it
        if not ok:                                        #If there is nothing in the image we skip the event  
            return False
        img,ROI=xtu.FindROI(img,ROI,self._roiwaistthres,self._roiexpand)                  #Crop the image, the ROI struct is changed. It also add an extra dimension to the image so the array can store multiple images corresponding to different bunches
        img=xtu.SplitImage(img,self._nb)
        imageStats=xtu.ProcessXTCAVImage(img,ROI)          #Obtain the different properties and profiles from the trace
                           
        PU=xtu.CalculatePhysicalUnits(ROI,[imageStats[0]['xCOM'],imageStats[0]['yCOM']],shotToShot,self._globalCalibration) #Obtain the physical units for the axis x and y, in fs and MeV
        
        #If the step in time is negative, we mirror the x axis to make it ascending and consequently mirror the profiles     
        if PU['xfsPerPix']<0:
            PU['xfs']=PU['xfs'][::-1]
            for j in range(self._nb):
                imageStats[j]['xProfile']=imageStats[j]['xProfile'][::-1]
                imageStats[j]['yCOMslice']=imageStats[j]['yCOMslice'][::-1]
                imageStats[j]['yRMSslice']=imageStats[j]['yRMSslice'][::-1]
                
        #Save the results of the step 1
        
        self._eventresultsstep1={
            'PU':PU,
            'imageStats':imageStats,
            'shotToShot':shotToShot
            }
        
        self._currenteventprocessedstep1=True
        return True
    
    def ProcessShotStep2(self):
        if not self._currenteventprocessedstep1:
            if not self.ProcessShotStep1():
                return False
        
        #Using all the available data, perform the retrieval for that given shot        
        self._eventresultsstep2=xtu.ProcessLasingSingleShot(self._eventresultsstep1['PU'],self._eventresultsstep1['imageStats'],self._eventresultsstep1['shotToShot'],self._lasingoffreference.averagedProfiles) 
        self._currenteventprocessedstep2=True  
        return True                    
            
    def InterBunchPulseDelayBasedOnCurrent(self):                               
        if not self._currenteventprocessedstep1:
            if not self.ProcessShotStep1():
                return []
            
        if (self._nb<2):
            return [0]
        
        t=self._eventresultsstep1['PU']['xfs']   
          
        peakpos=np.zeros((self._nb), dtype=np.float64);
        for j in range(0,self._nb):
            ind=np.mean(np.argpartition(-self._eventresultsstep1['imageStats'][j]['xProfile'],5)[0:5]) #Find the position of the 5 highest values
            peakpos[j]=t[ind]
            #peakpos[j]=t[np.argmax(self._eventresultsstep1['imageStats'][j]['xProfile'])]
                
        peakpos=peakpos-peakpos[0]
                    
        return peakpos
        
    def XRayPower(self):       
        
        t=[]
        power=[]
            
        if not self._currenteventprocessedstep2:
            if not self.ProcessShotStep2():
                return t,power
        
                        
        mastert=self._eventresultsstep2['t']
        t=np.zeros((self._nb,len(mastert)), dtype=np.float64);
        for j in range(0,self._nb):
            t[j,:]=mastert+self._eventresultsstep2['bunchdelay'][j]

        return t,(self._eventresultsstep2['powerERMS']+self._eventresultsstep2['powerECOM'])/2         
        
    def XRayPowerRMSBased(self):   
        
        t=[]
        power=[]
            
        if not self._currenteventprocessedstep2:
            if not self.ProcessShotStep2():
                return t,power
        
                        
        mastert=self._eventresultsstep2['t']
        t=np.zeros((self._nb,len(mastert)), dtype=np.float64);
        for j in range(0,self._nb):
            t[j,:]=mastert+self._eventresultsstep2['bunchdelay'][j]

        return t,self._eventresultsstep2['powerERMS']   
        
    def XRayPowerCOMBased(self):   
        
        t=[]
        power=[]
            
        if not self._currenteventprocessedstep2:
            if not self.ProcessShotStep2():
                return t,power
        
                        
        mastert=self._eventresultsstep2['t']
        t=np.zeros((self._nb,len(mastert)), dtype=np.float64);
        for j in range(0,self._nb):
            t[j,:]=mastert+self._eventresultsstep2['bunchdelay'][j]

        return t,self._eventresultsstep2['powerECOM']  
        
    def ReconstructionAgreement(self): 
    
        if not self._currenteventprocessedstep2:
            if not self.ProcessShotStep2():
                return float('nan')
                       
        return np.mean(self._eventresultsstep2['powerAgreement'])        
        
    @property
    def medianfilter(self):
        return self._medianfilter
    @medianfilter.setter
    def medianfilter(self, medianfilter):
        self._medianfilter = medianfilter
    @property
    def snrfilter(self):
        return self._snrfilter
    @snrfilter.setter
    def snrfilter(self, snrfilter):
        self._snrfilter = snrfilter
    @property
    def roiwaistthres(self):
        return self._roiwaistthres
    @roiwaistthres.setter
    def roiwaistthres(self, roiwaistthres):
        self._roiwaistthres = roiwaistthres
    @property
    def roiexpand(self):
        return self._groupsize
    @roiexpand.setter
    def roiexpand(self, roiexpand):
        self._roiexpand = roiexpand        
        
