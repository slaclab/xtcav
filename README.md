# XTCAV Analysis

XTCAV is a detector that is used to determine the laser-power vs. time of each LCLS shot.  Some detailed documentation from Tim Maxwell on this device is [here](https://confluence.slac.stanford.edu/display/PSDM/New+XTCAV+Documentation?preview=/181536250/181699034/xtcav-users-v0p4.pdf).  Alvaro Sanchez-Gonzalez authored the original psana-python code to do the rather complex analysis of images from the XTCAV camera to determine these quantities.  This current version has been updated to run more quickly and to fix some analysis errors.  

These scripts use some XTCAV data that was made public so they should be runnable by all users.  The scripts can be found in /reg/g/psdm/tutorials/examplePython/xtcav/ in the files xtcavDark.py, xtcavLasingOff.py, xtcavLasingOn.py.

## Quickstart

Two things must be done before XTCAV analysis will function: a "dark run" must be analyzed to get the pedestal values for cameras, and a "no lasing" run must be analyzed to generate sets of "no lasing" images (the latter is quite a complex process).  Note that for demonstration these first two scripts write constants to a "local" calibration directory called "calib".  For a real-experiment you won't need these lines because you will have permission to write to your official experiment calibration-constants directory.

Sample of dark run analysis:

```
# these two lines for example purposes only, to allow user to write
# calibration information to local directory called "calib"
# should be deleted for real analysis.

import psana
psana.setOption('psana.calib-dir','calib')
from xtcav.DarkBackgroundReference import *

DarkBackgroundReference(experiment='xpptut15', 
	run_number='300', 
	max_shots=1000,
	validity_range=(300,302))
```

Sample of "no-lasing" reference generating script:
```
# these two lines for example purposes only, to allow user to write
# calibration information to local directory called "calib"
# should be deleted for real analysis.
import psana
psana.setOption('psana.calib-dir','calib')
from xtcav.LasingOffReference import *
LasingOffReference(
	experiment='xpptut15',
	run_number='301',
	max_shots=1400,
	num_bunches=1,
	island_split_method = 'scipyLabel',
	validity_range=(300,302)) #only give first run number (i.e. (300,)) to have the validity range be open-ended ("end")

```


This script assumes that dark/lasing-off data has been analyzed (see above).  Unlike the previous two scripts it reads dark/lasing-off constants from the official calibration-directory. 

```
from psana import *
import matplotlib.pyplot as plt
from xtcav.LasingOnCharacterization import *
experiment = 'xpptut15'
run = '302'
mode = 'smd'
ds = psana.DataSource("exp=%s:run=%s:%s" % (experiment, run, mode))
ngood = 0

XTCAVRetrieval=LasingOnCharacterization() 

for evt in ds.events():
    if not XTCAVRetrieval.processEvent(evt):
        continue

    # time and power are lists, with each entry corresponding to
    # a bunch number.  The number of bunches is set by the GLOC.nb
    # parameter in the lasing-off analysis.  In general, one should
    # also cut on the "agreement" value, which measures the agreement
    # between the first and second moment analysis (larger is better).
    time, power = XTCAVRetrieval.xRayPower()  
    agreement = XTCAVRetrieval.reconstructionAgreement()

    ngood += 1
    print 'Agreement: %g%% ' % (agreement*100)
    plt.plot(time[0],power[0])
    plt.xlabel('Time (fs)')
    plt.ylabel('Lasing Power (GW)')
    plt.title('Agreement %4.2f'%agreement)
    plt.show()
    if ngood > 1: 
        break
    
```

* The LasingOnCharacerization module uses the [detector interface](https://confluence.slac.stanford.edu/pages/viewpage.action?pageId=205983617) to find the datasource being used. The program will fail if you try to process events without first setting a datasource.

* If you are analyzing an older experiment, you may find that psana does not support the 'smd' mode. Instead, use the 'idx' mode.


### Prerequisites

This code relies on the psana and PSCalib pacakages, which are automatically available on your SLAC UNIX account. Follow instructions [here](https://confluence.slac.stanford.edu/display/PSDM/psana+python+Setup) to make sure Beyond that, this package uses standard python packages such scipy, numpy, and cv2.


### Installing

If you would like to enhance or change the xtcav code, you can set up the repository locally using the usual 'git clone' method. However, you won't be able to test code on any of the psana data unless you're using one of the psana servers. For development, you should clone the repository into your UNIX account. Once you've ssh-ed into the psana server, in order to use the psana python package, you'll first have to run the command

```
source /reg/g/psdm/etc/psconda.sh
```
`cd` into the `xtcav` directory and run
```
python setup.py develop --user
```

This will allow the changes you make to persist without having to use a make file or reinstall the package. You can then test your changes and subsequently create a pull requests.

### Bugs & desired upgrades

The xtcav code is currently a work in progress. Some features missing include:

* Better electron bunch splitting capabilities
    * The current implementation should theoretically work for any number of bunches in the image. That being said, if bunches heavily intersect, then the code will only find one bunch. Recommendations for solving this include using the current 'connected components' finding method followed by a 'sparsest cut' method that would separate a single region into the two (or more) regions that minimizes the ratio (number of edges cut/min(size of regions)).

* Ability to have different numbers of 'clusters' for different electron bunches within same lasing off reference
    * Most of the groundwork for this has already been laid out. The current issue is saving the lasing off reference when the arrays are of variable length. This is an h5py problem. To fix this, the `FileInterface` file should be changed. Specifically, you'll need to use the special_dtype functionality in h5py to save lists of variable length arrays. The main problem is that it doesn't seem to work for varaible length arrays of arrays...
    Once this has been done, the line `num_groups = num_clusters` in `averageXTCAVProfileGroups` in `Utils.py` simply needs to be removed in order for this functionality to persist. 

* Lasing off profile clustering methods
    * Current clustering algorithms tested include Hierarchical (with cosine, l1, and euclidean distance metrics), KMeans, DBSCAN, and Birch. Analysis showed that Hierarchical with Euclidean affinity provided the best results. See child page of XTCAV confluence for specific results. The next step would be to compare these algorithms with the performance of a SVD composition method: <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4052871/>

* xtcavDisplay script
    * The scripts in `bin`, specifically `xtcavDisp`, are currently not very customizeable. Some desired features include providing flags to change the types of graphs shown, the size of graphs, etc. This would make it more useful for realtime anayses. 

* Add default values for xtcav global calibration values
    * The calibration values for the xtcav are currently populated from the psana datasource. It would be good to have some default values for variables such as umperpix, strstrength, etc. in case that information is missing or for running siumulations/experiments. Default values can be added to the `GlobalCalibration` namedtuple in Utils.py

* Variability in optimal number of clusters chosen
    * Because the reference sets used to calculate the gap statistic are generated randomly, there may be some variability in the number of clusters chosen. Therefore your results may be slightly different from run to run. If you’d like to avoid this, you can manually set the number of clusters chosen

