#Test GenerateDarkBackground
print "Testing GenerateDarkBackground"
from GenerateDarkBackground import *
GDB=GenerateDarkBackground();
GDB.experiment='amo86815'
GDB.runs='70'                       #The run number is irrelevant
GDB.maxshots=10
GDB.SetValidityRange(8000,8001)
GDB.Generate(savetofile=False);  #In the test we do everything but saving
print "Test GenerateDarkBackground successful"

#Test GenerateLasingOffReference
print "Testing GenerateLasingOffReference"
from GenerateLasingOffReference import *
GLOC=GenerateLasingOffReference();
GLOC.experiment='amo86815'
GLOC.runs='70'
GLOC.maxshots=11
GLOC.nb=2
GLOC.groupsize=5
GLOC.SetValidityRange(8000,8001)
splitmethods=['contourLabel','autothreshold']
for splitmethod in splitmethods:        
    GLOC.islandsplitmethod=splitmethod
    GLOC.Generate(savetofile=False);  #In the test we do everything but saving
    print "Test GenerateLasingOffReference successful for method: %s"%splitmethod


#Test ShotToShotCharacterization
print "Testing ShotToShotCharacterization"

import psana
from ShotToShotCharacterization import *

exp='amo86815'
runnum=70

dataSource=psana.DataSource("exp=%s:run=%d:idx"%(exp,runnum))

XTCAVRetrieval=ShotToShotCharacterization();
XTCAVRetrieval.SetDataSource(dataSource)
XTCAVRetrieval.nb=2
XTCAVRetrieval.snrfilter=3

run=dataSource.runs().next()
times=times = run.times()
evt=run.event(times[4])  

splitmethods=['contourLabel','autothreshold']


for splitmethod in splitmethods:        
    XTCAVRetrieval.islandsplitmethod=splitmethod
    XTCAVRetrieval.SetCurrentEvent(evt)              
           
    delaylist,ok=XTCAVRetrieval.PulseDelay(method='COM')
    delaylist,ok=XTCAVRetrieval.PulseDelay(method='RMS')
    delaylist,ok=XTCAVRetrieval.PulseDelay(method='RMSCOM')
    pulsefwhm,ok=XTCAVRetrieval.PulseFWHM(method='COM')   
    pulsefwhm,ok=XTCAVRetrieval.PulseFWHM(method='RMS')  
    pulsefwhm,ok=XTCAVRetrieval.PulseFWHM(method='RMSCOM')      

    delaylist,ok=XTCAVRetrieval.InterBunchPulseDelayBasedOnCurrent()   
    delaylist,ok=XTCAVRetrieval.InterBunchPulseDelayBasedOnCurrentMultiple(n=2,filterwith=20)
    delaylist,ok=XTCAVRetrieval.InterBunchPulseDelayBasedOnCurrentFourierFiltered()
    t,current,ok=XTCAVRetrieval.ElectronCurrentPerBunch()
    
    t,power,ok=XTCAVRetrieval.XRayPower(method='COM')
    t,power,ok=XTCAVRetrieval.XRayPower(method='RMS')
    t,power,ok=XTCAVRetrieval.XRayPower(method='RMSCOM')
    
    t,power,ok=XTCAVRetrieval.XRayPowerRMSBased()
    t,power,ok=XTCAVRetrieval.XRayPowerCOMBased()
    
    energy,ok=XTCAVRetrieval.XRayEnergyPerBunch(method='COM')
    energy,ok=XTCAVRetrieval.XRayEnergyPerBunch(method='RMS')
    energy,ok=XTCAVRetrieval.XRayEnergyPerBunch(method='RMSCOM')
    
    energy,ok=XTCAVRetrieval.XRayEnergyPerBunchCOMBased()
    energy,ok=XTCAVRetrieval.XRayEnergyPerBunchRMSBased()
    
    image,ok=XTCAVRetrieval.RawXTCAVImage()  
    image,ok=XTCAVRetrieval.ProcessedXTCAVImage()  
    roi,ok=XTCAVRetrieval.ProcessedXTCAVImageROI()  
    
    agreement,ok=XTCAVRetrieval.ReconstructionAgreement()
       
    fullresults,ok=XTCAVRetrieval.GetFullResults()  
    if ok:
        print ("All results obtained for method: %s"%splitmethod)             
    else:
        print ("Full results not calculated for method: %s"%splitmethod)       
        
print "Tests finished"