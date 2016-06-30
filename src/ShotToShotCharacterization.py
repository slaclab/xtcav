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
import UtilsPsana as xtup
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
        islandsplitmethod (str): island splitting algorithm. Set to 'scipylabel' or 'contourLabel'  The defaults parameter is then one used for the lasing off reference or 'scipylabel'.
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
        self._islandsplitmethod=''               #Method for island splitting
        self._islandsplitpar1=float('nan')
        self._islandsplitpar2=float('nan')
        self._currentevent=[]
        self._eventresultsstep1=[]
        self._eventresultsstep2=[]
        self._eventresultsstep3=[]
        self._env=[]
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
        self._envinfo=False      

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
            cp=CalibrationPaths(self._env,self._calpath)       
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
            cp=CalibrationPaths(self._env,self._calpath)     
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
        if not self._islandsplitmethod:
            self._islandsplitmethod=self._lasingoffreference.parameters.get('islandsplitmethod',self._lasingoffreference.parameters.get('islandSplitMethod','scipylabel'))#to account for the change name that the variable name suffered at some point, and make it compatible with older lasing off reference files            
        if math.isnan(self._islandsplitpar1):        
            self._islandsplitpar1=self._lasingoffreference.parameters.get('islandsplitpar1',self._lasingoffreference.parameters.get('par1',3.))          
        if math.isnan(self._islandsplitpar2):        
            self._islandsplitpar2=self._lasingoffreference.parameters.get('islandsplitpar2',self._lasingoffreference.parameters.get('par2',5.))               
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
        if not self._islandsplitmethod:
            self._islandsplitmethod='scipylabel'       
        if math.isnan(self._islandsplitpar1):        
            self._islandsplitpar1=3.0
        if math.isnan(self._islandsplitpar2):        
            self._islandsplitpar2=5.0
                           
    def SetCurrentEvent(self,evt):
        """
        Method that sets a psana event to be the current event. Only after setting an event it is possible to query for results such as X-Ray power, or pulse delay. On the other hand, the calculations to get the reconstruction will not be done until the information itself is requested, so the call to this method should be quite fast.

        Args:
            evt (psana event): relevant event to retrieve information form
            
        Returns:
            True: All the input form detectors necessary for a good reconstruction are present in the event. 
            False: The information from some detectors is missing for that event. It may still be possible to get information.
        """
        ebeam = evt.get(psana.Bld.BldDataEBeamV7,self.ebeam_data)   
        if not ebeam:
            ebeam = evt.get(psana.Bld.BldDataEBeamV6,self.ebeam_data)  
        if not ebeam:
            ebeam = evt.get(psana.Bld.BldDataEBeamV5,self.ebeam_data)  
        gasdetector=evt.get(psana.Bld.BldDataFEEGasDetEnergy,self.gasdetector_data) 
        if not gasdetector:
            gasdetector=evt.get(psana.Bld.BldDataFEEGasDetEnergyV1,self.gasdetector_data) 
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
        self.SetEnv(datasource.env())
    def SetEnv(self,env):
        """
        After creating an instance of the ShotToShotCharacterization class, it is necessary to pass the env object for that data that is being analysed.

        Args:
            env (Env object): Env object that is going to be used for the analysis.
            
        """
        self._env=env
        self._envinfo=False    
        
    def ProcessShotStep1(self):
        """
        Method that runs the first step of the reconstruction, which consists of getting statistics from the XTCAV trace. This method is called automatically and should not be called by the user unless he has a knowledge of the operation done by this class internally. 

        Returns: True if it was successful, False otherwise
        """

        if not self._currenteventavailable:
            return False
  
        #It is important that this is open first so the experiment name is set properly (important for loading references)   
        if not self._envinfo:
            self._experiment=self._env.experiment()
            epicsstore=self._env.epicsStore();
            self._globalCalibration,ok1=xtup.GetGlobalXTCAVCalibration(epicsstore)
            self._saturationValue = xtup.GetCameraSaturationValue(epicsstore)
            self._roixtcav,ok2=xtup.GetXTCAVImageROI(epicsstore) 
            if ok1 and ok2: #If the information is not good, we try next event
                self._envinfo=True
            else:
                return False

        #It is important that the lasing off reference is open first, because it may reset the lasing off reference that needs to be loaded        
        if not self._loadedlasingoffreference:
            self.LoadLasingOffReference()
        
        if not self._loadeddarkreference:
            self.LoadDarkReference()

        if np.max(self._rawimage)>=self._saturationValue : #Detection if the image is saturated, we skip if it is
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
        if ROI['xN']<3 or ROI['yN']<3:
            print 'ROI too small',ROI['xN'],ROI['yN']
            return False
        img=xtu.SplitImage(img,self._nb, self._islandsplitmethod,self._islandsplitpar1,self._islandsplitpar2)

         

        imageStats=xtu.ProcessXTCAVImage(img,ROI)          #Obtain the different properties and profiles from the trace        
        
        #Save the results of the step 1
        
        self._eventresultsstep1={
            'processedImage':img,
            'NB':img.shape[0],
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

        shotToShot,ok = xtup.ShotToShotParameters(self._ebeam,self._gasdetector) #Obtain the shot to shot parameters necessary for the retrieval of the x and y axis in time and energy units
        if not ok: #If the information is not good, we skip the event
            return False
                   
        imageStats=self._eventresultsstep1['imageStats'];
        ROI=self._eventresultsstep1['ROI']
                  
        PU, ok=xtu.CalculatePhysicalUnits(ROI,[imageStats[0]['xCOM'],imageStats[0]['yCOM']],shotToShot,self._globalCalibration) #Obtain the physical units for the axis x and y, in fs and MeV
        if not ok: #If the information is not good, we skip the event
            return False

        #If the step in time is negative, we mirror the x axis to make it ascending and consequently mirror the profiles     
        if PU['xfsPerPix']<0:
            PU['xfs']=PU['xfs'][::-1]
            for j in range(self._eventresultsstep1['NB']):
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

        #If the nubmer of bunches in the reference is not equal to the number of found bunches, we cannot reconstruct
        if self._eventresultsstep1['NB']!=self._nb:
            return False
        
        #Using all the available data, perform the retrieval for that given shot        
        self._eventresultsstep3=xtu.ProcessLasingSingleShot(self._eventresultsstep2['PU'],self._eventresultsstep2['imageStats'],self._eventresultsstep2['shotToShot'],self._lasingoffreference.averagedProfiles) 
        self._currenteventprocessedstep3=True  
        return True            
        
    def GetPhysicalUnitsResults(self):
        """
        Method which returns a dictionary based list with the physical units for the cropped image

        Returns: 
            out1: List with the results
                'yMeVPerPix':         Number of MeV per pixel for the vertical axis of the image
                'xfsPerPix':          Number of fs per pixel for the horizontal axis of the image
                'xfs':                Horizontal axis of the image in fs
                'yMeV':               Vertical axis of the image in MeV
            out2: True if the retrieval was successful, False otherwise. 
        """
    
        if not self._currenteventprocessedstep2:
            if not self.ProcessShotStep2():
                return None,False
        
        return self._eventresultsstep2['PU'],True                
        
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
            
    def PulseDelay(self,method='RMSCOM'):    
        """
        Method which returns the time of lasing for each bunch based on the x-ray reconstruction. They delays are referred to the center of mass of the total current. The order of the delays goes from higher to lower energy electron bunches.
        Args:
            method (str): method to use to obtain the power profile. 'RMS', 'COM' or 'RMSCOM' (Average of both)
        Returns: 
            out1: List of the delays for each bunch.
            out2: True if the retrieval was successful, False otherwise. 
        """
        if not self._currenteventprocessedstep3:
            if not self.ProcessShotStep3():
                return [],False
            
        if (self._eventresultsstep1['NB']<1):
            return np.zeros((self._eventresultsstep1['NB']), dtype=np.float64),True
        
                  
        peakpos=np.zeros((self._eventresultsstep1['NB']), dtype=np.float64);
        for j in range(0,self._eventresultsstep1['NB']):
            t=self._eventresultsstep3['t']+self._eventresultsstep3['bunchdelay'][j]
            if method=='RMS':
                power=self._eventresultsstep3['powerERMS'][j]
            elif method=='COM':
                power=self._eventresultsstep3['powerECOM'][j]
            elif method=='RMSCOM':
                power=(self._eventresultsstep3['powerECOM'][j]+self._eventresultsstep3['powerERMS'][j])/2
            else:
                return [],False        
            #quadratic fit around 5 pixels method
            central=np.argmax(power)
            try:
                fit=np.polyfit(t[central-2:central+3],power[central-2:central+3],2)
                peakpos[j]=-fit[1]/(2*fit[0])
            except:
                return [],False  
            
        return peakpos,True
            
    def PulseFWHM(self,method='RMSCOM'):    
        """
        Method which returns the FWHM of the pulse generated by each bunch in fs. It uses the power profile. The order of the widths goes from higher to lower energy electron bunches.
        Args:
            method (str): method to use to obtain the power profile. 'RMS', 'COM' or 'RMSCOM' (Average of both)
        Returns: 
            out1: List of the full widths half maximum for each bunch.
            out2: True if the retrieval was successful, False otherwise. 
        """
        if not self._currenteventprocessedstep3:
            if not self.ProcessShotStep3():
                return [],False
            
        if (self._eventresultsstep1['NB']<1):
            return np.zeros((self._eventresultsstep1['NB']), dtype=np.float64),True
        
                  
        peakwidth=np.zeros((self._eventresultsstep1['NB']), dtype=np.float64);
        for j in range(0,self._eventresultsstep1['NB']):
            t=self._eventresultsstep3['t']+self._eventresultsstep3['bunchdelay'][j]
            if method=='RMS':
                power=self._eventresultsstep3['powerERMS'][j]
            elif method=='COM':
                power=self._eventresultsstep3['powerECOM'][j]
            elif method=='RMSCOM':
                power=(self._eventresultsstep3['powerECOM'][j]+self._eventresultsstep3['powerERMS'][j])/2
            else:
                return [],False        
            #quadratic fit around 5 pixels method
            threshold=np.max(power)/2
            abovethrestimes=t[power>=threshold]
            dt=t[1]-t[0]
            peakwidth[j]=abovethrestimes[-1]-abovethrestimes[0]+dt
            
        return peakwidth,True
      
    def InterBunchPulseDelayBasedOnCurrent(self):    
        """
        Method which returns the time of lasing for each bunch based on the peak electron current on each bunch. A lasing off reference is not necessary for this retrieval. The delays are referred to the center of mass of the total current. The order of the delays goes from higher to lower energy electron bunches.

        Returns: 
            out1: List with the delay for each bunch.
            out2: True if the retrieval was successful, False otherwise. 
        """
        if not self._currenteventprocessedstep2:
            if not self.ProcessShotStep2():
                return [],False
            
        if (self._eventresultsstep1['NB']<1):
            return np.zeros((self._eventresultsstep1['NB']), dtype=np.float64),True
        
        t=self._eventresultsstep2['PU']['xfs']   
          
        peakpos=np.zeros((self._eventresultsstep1['NB']), dtype=np.float64);
        for j in range(0,self._eventresultsstep1['NB']):
            #highest value method
            #peakpos[j]=t[np.argmax(self._eventresultsstep1['imageStats'][j]['xProfile'])]
            
            #five highest values method
            #ind=np.mean(np.argpartition(-self._eventresultsstep2['imageStats'][j]['xProfile'],5)[0:5]) #Find the position of the 5 highest values
            #peakpos[j]=t[ind]
            
            #quadratic fit around 5 pixels method
            central=np.argmax(self._eventresultsstep1['imageStats'][j]['xProfile'])
            try:
                fit=np.polyfit(t[central-2:central+3],self._eventresultsstep1['imageStats'][j]['xProfile'][central-2:central+3],2)
                peakpos[j]=-fit[1]/(2*fit[0])
            except:
                return [],False  
            
        return peakpos,True
        
    def InterBunchPulseDelayBasedOnCurrentMultiple(self,n=1,filterwith=7):    
        """
        Method which returns multiple possible times of lasing for each bunch based on the peak electron current on each bunch. A lasing off reference is not necessary for this retrieval. The delays are referred to the center of mass of the total current. The order of the delays goes from higher to lower energy electron bunches. Then within each bunch the "n" delays are orderer from highest peak current yo lowest peak current.
        Args:
            n (int): number of possible times of lasing (peaks in the electron current) to find per bunch
            filterwith (float): Witdh of the peak that is removed before searching for the next peak in the same bunch
        Returns: 
            out1: List with a list of "n" delays for each bunch.
            out2: True if the retrieval was successful, False otherwise. 
        """
        if not self._currenteventprocessedstep2:
            if not self.ProcessShotStep2():
                return [],False
            
        if (self._eventresultsstep1['NB']<1):
            return np.zeros((self._eventresultsstep1['NB']), dtype=np.float64),True
        
        t=self._eventresultsstep2['PU']['xfs']   
          
        peakpos=np.zeros((self._eventresultsstep1['NB'],n), dtype=np.float64);
           
        for j in range(0,self._eventresultsstep1['NB']):
            profile=self._eventresultsstep1['imageStats'][j]['xProfile'].copy()
            for k in range(n):
                #highest value method
                #peakpos[j]=t[np.argmax(self._eventresultsstep1['imageStats'][j]['xProfile'])]
                
                #five highest values method
                #ind=np.mean(np.argpartition(-self._eventresultsstep2['imageStats'][j]['xProfile'],5)[0:5]) #Find the position of the 5 highest values
                #peakpos[j]=t[ind]
                
                #quadratic fit around 5 pixels method
                central=np.argmax(profile)
                try:
                    fit=np.polyfit(t[central-2:central+3],profile[central-2:central+3],2)
                    peakpos[j,k]=-fit[1]/(2*fit[0])
                    filter=1-np.exp(-(t-peakpos[j,k])**2/(filterwith/(2*np.sqrt(np.log(2))))**2)
                    profile=profile*filter                   
                except:
                    peakpos[j,k]=np.nan
                    if k==0:
                        return [],False
                
        return peakpos,True
        
    def InterBunchPulseDelayBasedOnCurrentFourierFiltered(self,targetwidthfs=20,thresholdfactor=0):    
        """
        Method which returns the time delay between the x-rays generated from different bunches based on the peak electron current on each bunch. A lasing off reference is not necessary for this retrieval. The delays are referred to the center of mass of the total current. The order of the delays goes from higher to lower energy electron bunches. This method includes a Fourier filter that applies a low pass filter to amplify the feature identified as the lasing part of the bunch, and ignore other peaks that may be higher in amplitude but also higher in width. It is possible to threshold the signal before calculating the Fourier transform to automatically discard peaks that may be sharp, but too low in amplitude to be the right peaks.
        Args:
            targetwidthfs (float): Witdh of the peak to be used for calculating delay
            thresholdfactor (float): Value between 0 and 1 that indicates which threshold factor to apply to filter the signal before calculating the fourier transform
        Returns: 
            out1: List with the delay for each bunch.
            out2: True if the retrieval was successful, False otherwise. 
        """
        if not self._currenteventprocessedstep2:
            if not self.ProcessShotStep2():
                return [],False
            
        if (self._eventresultsstep1['NB']<1):
            return np.zeros((self._eventresultsstep1['NB']), dtype=np.float64),True
        
        t=self._eventresultsstep2['PU']['xfs']   
        
        #Preparing the low pass filter
        N=len(t)
        dt=abs(self._eventresultsstep2['PU']['xfsPerPix'])
        if dt*N==0:
            return [],False
        df=1/(dt*N)
        
        f=np.array(range(0,N/2+1)+range(-N/2+1,0))*df
                           
        ffilter=(1-np.exp(-(f*targetwidthfs)**6))
          
        peakpos=np.zeros((self._eventresultsstep1['NB']), dtype=np.float64);
        for j in range(0,self._eventresultsstep1['NB']):
            #Getting the profile and the filtered version
            profile=self._eventresultsstep1['imageStats'][j]['xProfile']
            profilef=profile-np.max(profile)*thresholdfactor
            profilef[profilef<0]=0
            profilef=np.fft.ifft(np.fft.fft(profilef)*ffilter)
        
            #highest value method
            #peakpos[j]=t[np.argmax(profilef)]
            
            #five highest values method
            #ind=np.mean(np.argpartition(-profilef,5)[0:5]) #Find the position of the 5 highest values
            #peakpos[j]=t[ind]
            
            #quadratic fit around 5 pixels method and then fit to the original signal
            central=np.argmax(profilef)
            try:
                fit=np.polyfit(t[central-2:central+3],profile[central-2:central+3],2)
                peakpos[j]=-fit[1]/(2*fit[0])
            except:
                return [],False  
            
        return peakpos,True

    def QuadRefine(self,p):
        x1,x2,x3 = p + np.array([-1,0,1])
        y1,y2,y3 = self.wf[(p-self.rangelim[0]-1):(p-self.rangelim[0]+2)]
        d = (x1-x2)*(x1-x3)*(x2-x3)
        A = ( x3 * (y2-y1) + x2 * (y1-y3) + x1 * (y3-y2) ) / d
        B = ( x3**2.0 * (y1-y2) + x2**2.0 * (y3-y1) + x1**2.0 * (y2-y3) ) / d
        return -1*B / (2*A)

    def ElectronCurrentPerBunch(self):    
        """
        Method which returns the electron current per bunch. A lasing off reference is not necessary for this retrieval.

        Returns: 
            out1: time vectors in fs
            out2: electron currents in arbitrary units
            out3: True if the retrieval was successful, False otherwise
        """
        if not self._currenteventprocessedstep2:
            if not self.ProcessShotStep2():
                return [],[],False
        
        t=self._eventresultsstep2['PU']['xfs']   

        tout=np.zeros((self._eventresultsstep1['NB'],len(t)), dtype=np.float64);
        currents=np.zeros((self._eventresultsstep1['NB'],len(t)), dtype=np.float64);
        for j in range(0,self._eventresultsstep1['NB']):
            tout[j,:]=t
            currents[j,:]=self._eventresultsstep1['imageStats'][j]['xProfile']
                    
        return tout,currents,True
        
    def XRayPower(self,method='RMSCOM'):       
        """
        Method which returns the power profile for the X-Rays generated by each electron bunch. This is the averaged result from the RMS method and the COM method.

        Args:
            method (str): method to use to obtain the power profile. 'RMS', 'COM' or 'RMSCOM' (Average of both)
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

        if method=='RMS':
            power=self._eventresultsstep3['powerERMS']
        elif method=='COM':
            power=self._eventresultsstep3['powerECOM']
        elif method=='RMSCOM':
            power=(self._eventresultsstep3['powerECOM']+self._eventresultsstep3['powerERMS'])/2
        else:
            return t,[],False  
            
        return t,power,True         
        
    def XRayPowerRMSBased(self):   
        """
        Method which returns the power profile for the X-Rays generated by each electron bunch using the RMS method.

        Returns: 
            out1: time vectors in fs. 2D array where the first index refers to bunch number, and the second index to time.
            out2: power profiles in GW. 2D array where the first index refers to bunch number, and the second index to the power profile.
            out3: True if the retrieval was successful, False otherwise. 
        """
        
        return self.XRayPower(method='RMS') 

    def XRayPowerCOMBased(self):   
        """
        Method which returns the power profile for the X-Rays generated by each electron bunch using the COM method.

        Returns: 
            out1: time vectors in fs. 2D array where the first index refers to bunch number, and the second index to time.
            out2: power profiles in GW. 2D array where the first index refers to bunch number, and the second index to the power profile.
            out3: True if the retrieval was successful, False otherwise.
        """
        return self.XRayPower(method='COM') 
        
    def XRayEnergyPerBunch(self,method='RMSCOM'):   
        """
        Method which returns the total X-Ray energy generated per bunch. This is the averaged result from the RMS method and the COM method.
        Args:
            method (str): method to use to obtain the power profile. 'RMS', 'COM' or 'RMSCOM' (Average of both)
        Returns: 
            out1: List with the values of the energy for each bunch in J
            out2: True if the retrieval was successful, False otherwise.
        """
        energies=[]
            
        if not self._currenteventprocessedstep3:
            if not self.ProcessShotStep3():
                return energies,False
        
        
        if method=='RMS':
            energyperbunch=self._eventresultsstep3['lasingenergyperbunchERMS']
        elif method=='COM':
            energyperbunch=self._eventresultsstep3['lasingenergyperbunchECOM']
        elif method=='RMSCOM':
            energyperbunch=(self._eventresultsstep3['lasingenergyperbunchECOM']+self._eventresultsstep3['lasingenergyperbunchERMS'])/2
        else:
            return energies,False
       
        return energyperbunch,True  
        
    def XRayEnergyPerBunchCOMBased(self):   
        """
        Method which returns the total X-Ray energy generated per bunch based on the COM method.

        Returns: 
            out1: List with the values of the energy for each bunch in J
            out2: True if the retrieval was successful, False otherwise.
        """       
        return self.XRayEnergyPerBunch(method='COM')
    
    def XRayEnergyPerBunchRMSBased(self):   
        """
        Method which returns the total X-Ray energy generated per bunch based on the RMS method.

        Returns: 
            out1: List with the values of the energy for each bunch in J
            out2: True if the retrieval was successful, False otherwise.
        """
        return self.XRayEnergyPerBunch(method='RMS')      
        
        
        
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
        
    def ProcessedXTCAVImageROI(self):    
        """
        Method which returns the position of the processed XTCAV image within the whole CCD after background subtraction, noise removal, region of interest cropping and multiple bunch separation. This does not require a lasing off reference.

        Returns: 
            out1: Dictionary with the region of interest parameters.
            out2: True if the retrieval was successful, False otherwise.
        """     
        if not self._currenteventprocessedstep1:
            if not self.ProcessShotStep1():
                return [],False
            
        return self._eventresultsstep1['ROI'],True
        
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
    def nb(self):
        return self._nb
    @nb.setter
    def nb(self, nb):
        if not self._loadedlasingoffreference:
            self._nb = nb
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
    @property
    def islandsplitmethod(self):
        return self._islandsplitmethod
    @islandsplitmethod.setter
    def islandsplitmethod(self, islandsplitmethod):
        self._islandsplitmethod = islandsplitmethod 
    @property
    def islandsplitpar1(self):
        return self._islandsplitpar1
    @islandsplitpar1.setter
    def islandsplitpar1(self, islandsplitpar1):
        self._islandsplitpar1 = islandsplitpar1 
    @property
    def islandsplitpar2(self):
        return self._islandsplitpar2
    @islandsplitpar2.setter
    def islandsplitpar2(self, islandsplitpar2):
        self._islandsplitpar2 = islandsplitpar2
