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
import warnings
import Utils as xtu
from DarkBackground import *
from LasingOffReference import *
from CalibrationPaths import *


class ShotToShotCharacterization(object):

    def __init__(self):
        
        #Handle warnings
        warnings.filterwarnings('always',module='Utils',category=UserWarning)
        warnings.filterwarnings('ignore',module='Utils',category=RuntimeWarning, message="invalid value encountered in divide")
        
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
        self._eventresultsstep3=[]
        self._datasource=[]
        self._globalcalibration=[]
        self._roixtcav=[]
        self._calpath=''
        
        #Different flags
        self._loadeddarkreference=False
        self._loadedlasingoffreference=False
        self._currenteventavailable=False
        self._currenteventprocessedstep1=False   #Step one is pure processing of the trace, without lasing off referencing
        self._currenteventprocessedstep2=False   #Step two is calculating physical units
        self._currenteventprocessedstep3=False   #Step three is the processing of the profiles with respect to the reference
        self._datasourceinfo=False      

        #Camera and type for the xtcav images
        self.xtcav_camera = psana.Source('DetInfo(XrayTransportDiagnostic.0:Opal1000.0)')
        self.xtcav_type=psana.Camera.FrameV1
        self._rawimage=[]
        
        #Ebeam type: it should actually be the version 5 which is the one that contains xtcav stuff
        self.ebeam_data=psana.Source('BldInfo(EBeam)')
        self._ebeam=[]

        #Gas detectors for the pulse energies
        self.gasdetector_data=psana.Source('BldInfo(FEEGasDetEnergy)')
        self._gasdetector=[]
        
    def LoadDarkReference(self):
        if not self._darkreferencepath:
            cp=CalibrationPaths(self._datasource,self._calpath)       
            self._darkreferencepath=cp.findCalFileName('pedestals',self._currentevent.run())
            
        #If we could not find it, we just wont use it, and return False
        if not self._darkreferencepath:
            warnings.warn_explicit('Dark reference for run %d not found, image will not be background substracted' % self._currentevent.run(),UserWarning,'XTCAV',0)
            self._loadeddarkreference=False      
            return False
        
        self._darkreference=DarkBackground.Load(self._darkreferencepath)
        self._loadeddarkreference=True
        
        return True
        
    def LoadLasingOffReference(self):
        if not self._lasingoffreferencepath:
            cp=CalibrationPaths(self._datasource,self._calpath)     
            self._lasingoffreferencepath=cp.findCalFileName('lasingoffreference',self._currentevent.run())
            
        #If we could not find it, we load default parameters, and return False
        if not self._lasingoffreferencepath:
            warnings.warn_explicit('Lasing off reference for run %d not found, using set or default values for image processing' % self._currentevent.run(),UserWarning,'XTCAV',0)
            self.LoadDefaultProcessingParameters()            
            self._loadedlasingoffreference=False
            return False
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
            
        return True
            
    def LoadDefaultProcessingParameters(self):
        if math.isnan(self._nb):
            self._nb=1
        if math.isnan(self._medianfilter):
            self._medianfilter=3
        if math.isnan(self._snrfilter):
            self._snrfilter=10
        if math.isnan(self._roiwaistthres):
            self._roiwaistthres=0.2
        if math.isnan(self._roiexpand):
            self._roiexpand=2.5    
                           
    def SetCurrentEvent(self,evt):
        ebeam = evt.get(psana.Bld.BldDataEBeamV6,self.ebeam_data)   
        if not ebeam:
            ebeam = evt.get(psana.Bld.BldDataEBeamV5,self.ebeam_data)  
        gasdetector=evt.get(psana.Bld.BldDataFEEGasDetEnergy,self.gasdetector_data) 
        frame = evt.get(self.xtcav_type, self.xtcav_camera) 
        
        self._currenteventprocessedstep1=False
        self._currenteventprocessedstep2=False
        self._currenteventprocessedstep3=False
        self._eventresultsstep1=[]
        self._eventresultsstep2=[]
        self._eventresultsstep3=[]
        
        
        # If there is not frame, there is nothing we can do
        if (not frame):
            self._currenteventavailable=False     
            return False            
        
        self._rawimage=frame.data16().astype(np.float64)  
        self._currentevent=evt          
        self._ebeam=ebeam
        self._gasdetector=gasdetector        
        self._currenteventavailable=True
        # If gas detector or ebeam info is missing, we sill still may be able to do some stuff, but still return False
        if (ebeam and gasdetector):                               
            return True
        else:
            return False
        
    def SetDataSource(self,datasource):
        self._datasource=datasource
        self._datasourceinfo=False  

    def GetFullResults(self):
        if not self._currenteventprocessedstep3:
            if not self.ProcessShotStep3():
                return [],False
            
        return self._eventresultsstep3 ,True  
        
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
                                                            
        if np.max(self._rawimage)>=16383 : #Detection if the image is saturated, we skip if it is
            warnings.warn_explicit('Saturated Image',UserWarning,'XTCAV',0)
                                    
        #Subtract the dark background, taking into account properly possible different ROIs
        #Only if the reference is present
        if self._loadeddarkreference:        
            img,ROI=xtu.SubtractBackground(self._rawimage,self._roixtcav,self._darkreference.image,self._darkreference.ROI)  
        else:
            ROI=self._roixtcav
            img=self._rawimage
            
        img,ok=xtu.DenoiseImage(img,self._medianfilter,self._snrfilter)                    #Remove noise from the image and normalize it
        if not ok:                                        #If there is nothing in the image we skip the event  
            return False
        img,ROI=xtu.FindROI(img,ROI,self._roiwaistthres,self._roiexpand)                  #Crop the image, the ROI struct is changed. It also add an extra dimension to the image so the array can store multiple images corresponding to different bunches
        img=xtu.SplitImage(img,self._nb)
        imageStats=xtu.ProcessXTCAVImage(img,ROI)          #Obtain the different properties and profiles from the trace        
        
        #Save the results of the step 1
        
        self._eventresultsstep1={
            'processedImage':img,
            'ROI':ROI,
            'imageStats':imageStats,
            }
        
        self._currenteventprocessedstep1=True
                        
        return True
        
    def ProcessShotStep2(self):
        if not self._currenteventprocessedstep1:
            if not self.ProcessShotStep1():
                return False
        
        shotToShot,ok = xtu.ShotToShotParameters(self._ebeam,self._gasdetector) #Obtain the shot to shot parameters necessary for the retrieval of the x and y axis in time and energy units
        if not ok: #If the information is not good, we skip the event
            return False
                           
        imageStats=self._eventresultsstep1['imageStats'];
        ROI=self._eventresultsstep1['ROI']
                           
        PU=xtu.CalculatePhysicalUnits(ROI,[imageStats[0]['xCOM'],imageStats[0]['yCOM']],shotToShot,self._globalCalibration) #Obtain the physical units for the axis x and y, in fs and MeV
        
        #If the step in time is negative, we mirror the x axis to make it ascending and consequently mirror the profiles     
        if PU['xfsPerPix']<0:
            PU['xfs']=PU['xfs'][::-1]
            for j in range(self._nb):
                imageStats[j]['xProfile']=imageStats[j]['xProfile'][::-1]
                imageStats[j]['yCOMslice']=imageStats[j]['yCOMslice'][::-1]
                imageStats[j]['yRMSslice']=imageStats[j]['yRMSslice'][::-1]
                
        #Save the results of the step 2
        
        self._eventresultsstep2={
            'PU':PU,
            'imageStats':imageStats,
            'shotToShot':shotToShot
            }
        
        self._currenteventprocessedstep2=True
        return True       
    
    def ProcessShotStep3(self):
        if not self._currenteventprocessedstep2:
            if not self.ProcessShotStep2():
                return False
        
        #There is no possible step 3 if there is not lasing off reference
        if not self._loadedlasingoffreference:
            return False
        
        #Using all the available data, perform the retrieval for that given shot        
        self._eventresultsstep3=xtu.ProcessLasingSingleShot(self._eventresultsstep2['PU'],self._eventresultsstep2['imageStats'],self._eventresultsstep2['shotToShot'],self._lasingoffreference.averagedProfiles) 
        self._currenteventprocessedstep3=True  
        return True                    
            
    def InterBunchPulseDelayBasedOnCurrent(self):                               
        if not self._currenteventprocessedstep2:
            if not self.ProcessShotStep2():
                return [],False
            
        if (self._nb<2):
            return [0],True
        
        t=self._eventresultsstep2['PU']['xfs']   
          
        peakpos=np.zeros((self._nb), dtype=np.float64);
        for j in range(0,self._nb):
            ind=np.mean(np.argpartition(-self._eventresultsstep2['imageStats'][j]['xProfile'],5)[0:5]) #Find the position of the 5 highest values
            peakpos[j]=t[ind]
            #peakpos[j]=t[np.argmax(self._eventresultsstep1['imageStats'][j]['xProfile'])]
                
        peakpos=peakpos-peakpos[0]
                    
        return peakpos,True
        
    def XRayPower(self):       
        
        t=[]
        power=[]
            
        if not self._currenteventprocessedstep3:
            if not self.ProcessShotStep3():
                return t,power,False
        
                        
        mastert=self._eventresultsstep3['t']
        t=np.zeros((self._nb,len(mastert)), dtype=np.float64);
        for j in range(0,self._nb):
            t[j,:]=mastert+self._eventresultsstep3['bunchdelay'][j]

        return t,(self._eventresultsstep3['powerERMS']+self._eventresultsstep3['powerECOM'])/2,True         
        
    def XRayPowerRMSBased(self):   
        
        t=[]
        power=[]
            
        if not self._currenteventprocessedstep3:
            if not self.ProcessShotStep3():
                return t,power,False
        
                        
        mastert=self._eventresultsstep3['t']
        t=np.zeros((self._nb,len(mastert)), dtype=np.float64);
        for j in range(0,self._nb):
            t[j,:]=mastert+self._eventresultsstep3['bunchdelay'][j]

        return t,self._eventresultsstep3['powerERMS'],True   
        

    def XRayPowerCOMBased(self):   
        
        t=[]
        power=[]
            
        if not self._currenteventprocessedstep3:
            if not self.ProcessShotStep3():
                return t,power,False
        
                        
        mastert=self._eventresultsstep3['t']
        t=np.zeros((self._nb,len(mastert)), dtype=np.float64);
        for j in range(0,self._nb):
            t[j,:]=mastert+self._eventresultsstep3['bunchdelay'][j]

        return t,self._eventresultsstep3['powerECOM'],True  
        
    def XRayEnergyPerBunch(self):   
        
        energies=[]
            
        if not self._currenteventprocessedstep3:
            if not self.ProcessShotStep3():
                return energies,False
       
        return (self._eventresultsstep3['lasingenergyperbunchECOM']+self._eventresultsstep3['lasingenergyperbunchERMS'])/2,True  
        
    def XRayEnergyPerBunchCOMBased(self):   
        
        energies=[]
            
        if not self._currenteventprocessedstep3:
            if not self.ProcessShotStep3():
                return energies,False
       
        return self._eventresultsstep3['lasingenergyperbunchECOM'],True  
    
    def XRayEnergyPerBunchRMSBased(self):   
        
        energies=[]
            
        if not self._currenteventprocessedstep3:
            if not self.ProcessShotStep3():
                return energies,False
       
        return self._eventresultsstep3['lasingenergyperbunchERMS'],True          
        
        
        
    def RawXTCAVImage(self):                  
        if not self._currenteventavailable:
            return [],False
            
        return self._rawimage,True
        
    def ProcessedXTCAVImage(self):                  
        if not self._currenteventprocessedstep1:
            if not self.ProcessShotStep1():
                return [],False
            
        return self._eventresultsstep1['processedImage'],True
        
    def ReconstructionAgreement(self): 
    
        if not self._currenteventprocessedstep3:
            if not self.ProcessShotStep3():
                return float('nan'),False
                       
        return np.mean(self._eventresultsstep3['powerAgreement'])  ,True      
        
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
    @property
    def calibrationpath(self):
        return self._calpath
    @calibrationpath.setter
    def calibrationpath(self, calpath):
        self._calpath = calpath       
        
