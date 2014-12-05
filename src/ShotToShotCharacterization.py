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

    """
    Class that can be used to reconstruct the full X-Ray power time profile for single or multiple bunches, relying on the presence of a dark background reference, and a lasing off reference. (See GenerateDarkBackground and Generate LasingOffReference for more information)
    Attributes:
        calibrationpath (str): Custom calibration directory in case the default is not intended to be used.
        medianfilter (int): Number of neighbours for median filter (If not set, the value that was used for the lasing off reference will be used).
        snrfilter (float): Number of sigmas for the noise threshold (If not set, the value that was used for the lasing off reference will be used).
        roiwaistthres (float): ratio with respect to the maximum to decide on the waist of the XTCAV trace (If not set, the value that was used for the lasing off reference will be used).
        roiexpand (float): number of waists that the region of interest around will span around the center of the trace (If not set, the value that was used for the lasing off reference will be used).
    """

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
        """
        Method that loads the dark reference. This method is called automatically and should not be called by the user unless he has a knowledge of the operation done by this class internally.
        
        Returns: True if successful, False otherwise        
        """
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
        """
        Method that loads the lasing off reference. This method is called automatically and should not be called by the user unless he has a knowledge of the operation done by this class internally.
        
        Returns: True if successful, False otherwise        
        """
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
        """
        Method that sets some standard processing parameters in case they have not been explicitly set by the user and could not been retrieved from the lasing off reference. This method is called automatically and should not be called by the user unless he has a knowledge of the operation done by this class internally.             
        """
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
        """
        Method that sets a psana event to be the current event. Only after setting an event it is possible to query for results such as X-Ray power, or pulse delay. On the other hand, the calculations to get the reconstruction will not be done until the information itself is requested, so the call to this method should be quite fast.

        Args:
            evt (psana event): relevant event to retrieve information form
            
        Returns:
            True: All the input form detectors necessary for a good reconstruction are present in the event. 
            False: The information from some detectors is missing for that event. It may still be possible to get information.
        """
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
        """
        After creating an instance of the ShotToShotCharacterization class, it is necessary to pass the dataSource object for that data that is being analysed.

        Args:
            datasource (Datasource object): DataSource object that is going to be used for the analysis.
            
        """
        self._datasource=datasource
        self._datasourceinfo=False    
        
    def ProcessShotStep1(self):
        """
        Method that runs the first step of the reconstruction, which consists of getting statistics from the XTCAV trace. This method is called automatically and should not be called by the user unless he has a knowledge of the operation done by this class internally. 

        Returns: True if it was successful, False otherwise
        """
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
        """
        Method that runs the second step of the reconstruction, which consists of converting from pixel units into time and energy units for the trace. This method is called automatically and should not be called by the user unless he has a knowledge of the operation done by this class internally. 

        Returns: True if it was successful, False otherwise
        """
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
        """
        Method that runs the third step of the reconstruction, which consists of comparing the profiles to the reference profiles to obtain the X-Ray power. This method is called automatically and should not be called by the user unless he has a knowledge of the operation done by this class internally. 

        Returns: True if it was successful, False otherwise
        """
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

    def GetFullResults(self):
        """
        Method which returns a dictionary based list with the full results of the characterization

        Returns: 
            out1: List with the results
                 't':                           Master time vector in fs
                'powerECOM':                    Retrieved power in GW based on ECOM
                'powerERMS':                    Retrieved power in GW based on ERMS
                'powerAgreement':               Agreement between the two intensities
                'bunchdelay':                   Delay from each bunch with respect to the first one in fs
                'bunchdelaychange':             Difference between the delay from each bunch with respect to the first one in fs and the same form the non lasing reference
                'xrayenergy':                   Total x-ray energy from the gas detector in J
                'lasingenergyperbunchECOM':     Energy of the XRays generated from each bunch for the center of mass approach in J
                'lasingenergyperbunchERMS':     Energy of the XRays generated from each bunch for the dispersion approach in J
                'bunchenergydiff':              Distance in energy for each bunch with respect to the first one in MeV
                'bunchenergydiffchange':        Comparison of that distance with respect to the no lasing
                'lasingECurrent':               Electron current for the lasing trace (In #electrons/s)
                'nolasingECurrent':             Electron current for the no lasing trace (In #electrons/s)
                'lasingECOM':                   Lasing energy center of masses for each time in MeV
                'nolasingECOM':                 No lasing energy center of masses for each time in MeV
                'lasingERMS':                   Lasing energy dispersion for each time in MeV
                'nolasingERMS':                 No lasing energy dispersion for each time in MeV
                'NB':                           Number of bunches
            out2: True if the retrieval was successful, False otherwise. 
        """
        if not self._currenteventprocessedstep3:
            if not self.ProcessShotStep3():
                return [],False
            
        return self._eventresultsstep3 ,True        
            
    def InterBunchPulseDelayBasedOnCurrent(self):    
        """
        Method which returns the time delay between the x-rays generated from different bunches with respect to the first one based on the peak electron current on each bunch. A lasing off reference is not necessary for this retrieval.

        Returns: 
            out1: List of the delays. Since it is referenced to the first bunch, the first value is always 0.
            out2: True if the retrieval was successful, False otherwise. 
        """
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
        """
        Method which returns the power profile for the X-Rays generated by each electron bunch. This is the averaged result from the RMS method and the COM method.

        Returns: 
            out1: time vectors in fs. 2D array where the first index refers to bunch number, and the second index to time.
            out2: power profiles in GW. 2D array where the first index refers to bunch number, and the second index to the power profile.
            out3: True if the retrieval was successful, False otherwise. 
        """
        
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
        """
        Method which returns the power profile for the X-Rays generated by each electron bunch using the RMS method.

        Returns: 
            out1: time vectors in fs. 2D array where the first index refers to bunch number, and the second index to time.
            out2: power profiles in GW. 2D array where the first index refers to bunch number, and the second index to the power profile.
            out3: True if the retrieval was successful, False otherwise. 
        """
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
        """
        Method which returns the power profile for the X-Rays generated by each electron bunch using the COM method.

        Returns: 
            out1: time vectors in fs. 2D array where the first index refers to bunch number, and the second index to time.
            out2: power profiles in GW. 2D array where the first index refers to bunch number, and the second index to the power profile.
            out3: True if the retrieval was successful, False otherwise.
        """
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
        """
        Method which returns the total X-Ray energy generated per bunch. This is the averaged result from the RMS method and the COM method.

        Returns: 
            out1: List with the values of the energy for each bunch in J
            out2: True if the retrieval was successful, False otherwise.
        """
        energies=[]
            
        if not self._currenteventprocessedstep3:
            if not self.ProcessShotStep3():
                return energies,False
       
        return (self._eventresultsstep3['lasingenergyperbunchECOM']+self._eventresultsstep3['lasingenergyperbunchERMS'])/2,True  
        
    def XRayEnergyPerBunchCOMBased(self):   
        """
        Method which returns the total X-Ray energy generated per bunch based on the COM method.

        Returns: 
            out1: List with the values of the energy for each bunch in J
            out2: True if the retrieval was successful, False otherwise.
        """
        energies=[]
            
        if not self._currenteventprocessedstep3:
            if not self.ProcessShotStep3():
                return energies,False
       
        return self._eventresultsstep3['lasingenergyperbunchECOM'],True  
    
    def XRayEnergyPerBunchRMSBased(self):   
        """
        Method which returns the total X-Ray energy generated per bunch based on the RMS method.

        Returns: 
            out1: List with the values of the energy for each bunch in J
            out2: True if the retrieval was successful, False otherwise.
        """
        energies=[]
            
        if not self._currenteventprocessedstep3:
            if not self.ProcessShotStep3():
                return energies,False
       
        return self._eventresultsstep3['lasingenergyperbunchERMS'],True          
        
        
        
    def RawXTCAVImage(self):     
        """
        Method which returns the raw XTCAV image. This does not require of references at all.

        Returns: 
            out1: 2D array with the image
            out2: True if the retrieval was successful, False otherwise.
        """    
        if not self._currenteventavailable:
            return [],False
            
        return self._rawimage,True
        
    def ProcessedXTCAVImage(self):    
        """
        Method which returns the processed XTCAV image after background subtraction, noise removal, region of interest cropping and multiple bunch separation. This does not require a lasing off reference.

        Returns: 
            out1: 3D array where the first index is bunch number, and the other two are the image.
            out2: True if the retrieval was successful, False otherwise.
        """     
        if not self._currenteventprocessedstep1:
            if not self.ProcessShotStep1():
                return [],False
            
        return self._eventresultsstep1['processedImage'],True
        
    def ReconstructionAgreement(self): 
        """
        Value for the agreement of the reconstruction using the RMS method and using the COM method. It consists of a value ranging from -1 to 1.

        Returns: 
            out1: value for the agreement.
            out2: True if the retrieval was successful, False otherwise.
        """
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
        return self._roiexpand
    @roiexpand.setter
    def roiexpand(self, roiexpand):
        self._roiexpand = roiexpand
    @property
    def calibrationpath(self):
        return self._calpath
    @calibrationpath.setter
    def calibrationpath(self, calpath):
        self._calpath = calpath       
        
