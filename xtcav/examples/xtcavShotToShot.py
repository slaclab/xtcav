import psana
from xtcav.ShotToShotCharacterization import *
import numpy as np

runs = '137'
experiment='amox23616'
maxshots = 5

XTCAVRetrieval=ShotToShotCharacterization() 

#Process individual events
data_source = psana.DataSource("exp=%s:run=%s:idx" % (experiment, runs))

XTCAVRetrieval.SetDataSource(data_source)
for run in data_source.runs():
    n_r=0  #Counter for the total number of xtcav images processed within the run       
    times = run.times()
    for t in times:
        evt = run.event(t)
        if not XTCAVRetrieval.processEvent(evt):
            continue

        t, power = XTCAVRetrieval.XRayPower()  
        agreement = XTCAVRetrieval.ReconstructionAgreement()
        pulse = XTCAVRetrieval.PulseDelay()
        print 'Agreement: %g%% Maximum power: %g GW Pulse Delay: %g ' %(agreement*100,np.amax(power), pulse[0])
        # image = XTCAVRetrieval.ProcessedXTCAVImage()
        # np.save("testing", image)
        n_r += 1            

        if n_r>=maxshots: 
            break
