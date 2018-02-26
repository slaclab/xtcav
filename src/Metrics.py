import collections
import Constants
import numpy as np

"""
    Wrappers for metrics for XTCAV analysis. 
    Most perform checks to make sure information is valid. Using 1/0 to denote
    validity since h5py does not support boolean values (at least that's how it appears..)
    """
def namedtuple(typename, field_names, default_values=()):
    T = collections.namedtuple(typename, field_names)
    T.__new__.__defaults__ = (None,) * len(T._fields)

    if isinstance(default_values, collections.Mapping):
        prototype = T(**default_values)
    else:
        prototype = T(*default_values)

    T.__new__.__defaults__ = tuple(prototype)
    return T


ROIMetrics = namedtuple('ROIMetrics',
    ['xN', #Size of the image in X   
    'x0',  #Position of the first pixel in x
    'yN',  #Size of the image in Y 
    'y0',  #Position of the first pixel in y
    'x',   #X vector
    'y',   #Y vector
    'valid'], 
    {'valid': 0,
     'xN': 1024,                      
     'x0': 0, 
     'yN': 1024, 
     'y0': 0,
     'x': np.arange(0, 1024),
     'y': np.arange(0, 1024)})


GlobalCalibration = namedtuple('GlobalCalibration', 
    ['umperpix', #Pixel size of the XTCAV camera
    'strstrength', #Strength parameter
    'rfampcalib', #Calibration of the RF amplitude
    'rfphasecalib', #Calibration of the RF phase
    'dumpe',        #Beam energy: dump config
    'dumpdisp'])

      
ShotToShotParameters = namedtuple('ShotToShotParameters',
    ['ebeamcharge',  #ebeamcharge
    'dumpecharge',  #dumpecharge in C
    'xtcavrfamp',   #RF amplitude
    'xtcavrfphase', #RF phase
    'xrayenergy',   #Xrays energy in J
    'unixtime',
    'fiducial',
    'valid'],
    {'ebeamcharge': Constants.E_BEAM_CHARGE,
    'dumpecharge': Constants.DUMP_E_CHARGE,
    'xtcavrfphase': Constants.XTCAV_RFPHASE,
    'xtcavrfamp': Constants.XTCAV_RFAMP,
    'valid': 1}
    )


ImageStatistics = namedtuple('ImageStatistics', 
    ['imfrac',
    'xProfile',
    'yProfile',
    'xCOM',
    'yCOM',
    'xRMS',
    'yRMS',
    'xFWHM',
    'yFWHM',
    'yCOMslice',
    'yRMSslice'])


PhysicalUnits = namedtuple('PhysicalUnits', 
    ['xfs',
    'yMeV',
    'xfsPerPix',
    'yMeVPerPix',
    'valid'])

ImageProfile = namedtuple('ImageProfile', 
    ['image_stats',
    'roi',
    'shot_to_shot',
    'physical_units'])


LasingOffParameters = namedtuple('LasingOffParameters', 
    ['experiment', 
    'maxshots', 
    'run', 
    'start',
    'validityrange', 
    'darkreferencepath', 
    'num_bunches', 
    'groupsize', 
    'medianfilter', 
    'snrfilter', 
    'roiwaistthres', 
    'roiexpand', 
    'islandsplitmethod',
    'islandsplitpar1', 
    'islandsplitpar2', 
    'calpath', 
    'version'])


DarkBackgroundParameters = namedtuple('DarkBackgroundParameters', 
    ['experiment', 
    'maxshots', 
    'run', 
    'validityrange', 
    'calibrationpath'])

AveragedProfiles = namedtuple('AveragedProfiles',
    ['t',                         #Master time in fs
    'eCurrent',                   #Electron current in (#electrons/s)
    'eCOMslice',                  #Energy center of masses for each time in MeV
    'eRMSslice',                  #Energy dispersion for each time in MeV
    'distT',                      #Distance in time of the center of masses with respect to the center of the first bunch in fs
    'distE',                      #Distance in energy of the center of masses with respect to the center of the first bunch in MeV
    'tRMS',                       #Total dispersion in time in fs
    'eRMS',                       #Total dispersion in energy in MeV
    'num_bunches',                         #Number of bunches
    'num_groups',                          #Number of profiles
    'eventTime',                  #Unix times used for jumping to events
    'eventFid'])                  #Fiducial values used for jumping to events


