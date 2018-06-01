import psana
from xtcav.LasingOnCharacterization import *
import numpy as np

runs = '574'
experiment='diamcc14'
maxshots = 50

XTCAVRetrieval=LasingOnCharacterization() 

#Process individual events
# data_source = psana.DataSource("exp=%s:run=%s:idx" % (experiment, runs))
# n_r=0  #Counter for the total number of xtcav images processed within the run 
# run = data_source.runs().next()

# times = run.times()
# num_processed = 0 #Counter for the total number of xtcav images processed within the run  
# for t in times: 
#     evt = run.event(t)
#     if not XTCAVRetrieval.processEvent(evt):
#         continue

#     t, power = XTCAVRetrieval.xRayPower()  
#     agreement = XTCAVRetrieval.reconstructionAgreement()
#     pulse = XTCAVRetrieval.pulseDelay()
#     print 'Agreement: %g%% Maximum power: %g GW Pulse Delay: %g ' %(agreement*100,np.amax(power), pulse[0])
    
#     n_r += 1   

#     if n_r>=maxshots: 
#         break
data_source = psana.DataSource("exp=%s:run=%s:smd" % ('amox23616', '138'))
n_r=0  #Counter for the total number of xtcav images processed within the run 
arr=[]
arr2=[]
for evt in data_source.events():
    if not XTCAVRetrieval.processEvent(evt):
        continue

    t, power = XTCAVRetrieval.xRayPower(method='RMS')
    arr.append(power)  
    agreement = XTCAVRetrieval.reconstructionAgreement()
    pulse = XTCAVRetrieval.pulseDelay()
    print 'Agreement: %g%% Maximum power: %g GW Pulse Delay: %g ' %(agreement*100,np.amax(power), pulse[0])
    arr2.append(XTCAVRetrieval.processedXTCAVImage()[0])
    n_r += 1   

    if n_r>=maxshots: 
        break
np.save("images", arr2)
np.save("power", arr)
