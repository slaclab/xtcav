import psana
from xtcav.LasingOnCharacterization import *
import numpy as np

runs = '82'
experiment='amox23616'
maxshots = 50

XTCAVRetrieval=LasingOnCharacterization(lasingoffreferencepath='/reg/d/psdm/AMO/amox23616/calib/Xtcav::CalibV1/XrayTransportDiagnostic.0:Opal1000.0/lasingoffreference/79-end.data') 

#Process individual events
data_source = psana.DataSource("exp=%s:run=%s:idx" % (experiment, runs))
n_r=0  #Counter for the total number of xtcav images processed within the run 
arr = []
arr2 = []
run = data_source.runs().next()

times = run.times()
num_processed = 0 #Counter for the total number of xtcav images processed within the run  
arr = []
for t in times: 
    evt = run.event(t)
    if not XTCAVRetrieval.processEvent(evt):
        continue

    t, power = XTCAVRetrieval.xRayPower()  
    agreement = XTCAVRetrieval.reconstructionAgreement()
    pulse = XTCAVRetrieval.pulseDelay()
    print 'Agreement: %g%% Maximum power: %g GW Pulse Delay: %g ' %(agreement*100,np.amax(power), pulse[0])
    arr.append(XTCAVRetrieval.processedXTCAVImage()[0])
    arr2.append(XTCAVRetrieval.processedXTCAVImageProfile())
    
    n_r += 1   

    if n_r>=maxshots: 
        break

# n_r=0  #Counter for the total number of xtcav images processed within the run 
# arr = []
# arr2 = []
# for evt in data_source.events():
#     if not XTCAVRetrieval.processEvent(evt):
#         continue

#     t, power = XTCAVRetrieval.xRayPower()  
#     agreement = XTCAVRetrieval.reconstructionAgreement()
#     pulse = XTCAVRetrieval.pulseDelay()
#     print 'Agreement: %g%% Maximum power: %g GW Pulse Delay: %g ' %(agreement*100,np.amax(power), pulse[0])
#     arr.append(XTCAVRetrieval.processedXTCAVImage()[0])
#     arr2.append(XTCAVRetrieval.processedXTCAVImageProfile())
    
#     n_r += 1   

#     if n_r>=maxshots: 
#         break
np.save("lasingon_ex", arr)
np.save("lasingon_profile_ex", arr2)
