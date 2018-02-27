import psana
from xtcav.ShotToShotCharacterization import *

runs = '137'
experiment='amox23616'
maxshots = 5

XTCAVRetrieval=ShotToShotCharacterization(maxshots = maxshots,  #max number of valid shots to process
                                          experiment=experiment,
                                          runs = runs) 
# Process all runs at once
# processed_runs = XTCAVRetrieval.processRuns() 

#Process individual events
data_source  = psana.DataSource("exp=%s:run=%s:idx" % (experiment, runs))
for r, run in enumerate(data_source.runs()):
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

        n_r += 1            

        if n_r>=maxshots: 
            break
