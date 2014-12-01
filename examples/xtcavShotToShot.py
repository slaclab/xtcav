import psana
from xtcav.ShotToShotCharacterization import *

maxshots=5             #Maximum number of valid shots to process
experiment='amoc8114'  #Experiment label
runs='87'              #Runs

#Loading the dataset from the "dark" run, this way of working should be compatible with both xtc and hdf5 files
dataSource=psana.DataSource("exp=%s:run=%s:idx" % (experiment,runs))

#XTCAV Retrieval (setting the data source is useful to get information such as experiment name)
XTCAVRetrieval=ShotToShotCharacterization();
XTCAVRetrieval.SetDataSource(dataSource)

for r,run in enumerate(dataSource.runs()):
    n_r=0  #Counter for the total number of xtcav images processed within the run       
    times = run.times()
    for t in times:
        evt = run.event(t)

        if not XTCAVRetrieval.SetCurrentEvent(evt):
            continue

        t,power,ok1=XTCAVRetrieval.XRayPower()  
        agreement,ok2=XTCAVRetrieval.ReconstructionAgreement()      

        if (ok1 and ok2):
            print "%d/%d" % (n_r,maxshots) #Debugging purposes, will be removed
            print 'Agreement: %g %% Maximum power: %g GW' %(agreement*100,np.amax(power))

        n_r=n_r+1            

        if n_r>=maxshots: #After a certain number of shots we stop (Ideally this would be an argument, rather than a hardcoded value)
            break
