import psana
from xtcav.LasingOnCharacterization import *
import numpy as np

runs = '137'
experiment='amox23616'
maxshots = 50

XTCAVRetrieval=LasingOnCharacterization(lasingoffreferencepath="/reg/d/psdm/AMO/amox23616/calib/Xtcav::CalibV1/XrayTransportDiagnostic.0:Opal1000.0/lasingoffreference/131-end.data") 

#Process individual events
data_source = psana.DataSource("exp=%s:run=%s:smd" % (experiment, runs))

n_r=0  #Counter for the total number of xtcav images processed within the run 
for evt in data_source.events():
    if not XTCAVRetrieval.processEvent(evt):
        continue

    t, power = XTCAVRetrieval.xRayPower()  
    agreement = XTCAVRetrieval.reconstructionAgreement()
    pulse = XTCAVRetrieval.pulseDelay()
    print 'Agreement: %g%% Maximum power: %g GW Pulse Delay: %g ' %(agreement*100,np.amax(power), pulse[0])
    #image = XTCAVRetrieval.RawXTCAVImage()
    #np.save("test_img", image)
    n_r += 1   

    if n_r>=maxshots: 
        break
