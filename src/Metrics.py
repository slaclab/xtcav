import collections

"""
    Wrappers for metrics for XTCAV analysis. 
    Most perform checks to make sure information is valid. Using 1/0 to denote
    validity since h5py does not support boolean values (at least that's how it appears..)
    """
def namedtuple_with_defaults(typename, field_names, default_values=()):
    T = collections.namedtuple(typename, field_names)
    T.__new__.__defaults__ = (None,) * len(T._fields)
    if isinstance(default_values, collections.Mapping):
        prototype = T(**default_values)
    else:
        prototype = T(*default_values)
    T.__new__.__defaults__ = tuple(prototype)
    return T


ROIMetrics = namedtuple_with_defaults('ROIMetrics', 
    ['xN', #Size of the image in X   
    'x0',  #Position of the first pixel in x
    'yN',  #Size of the image in Y 
    'y0',  #Position of the first pixel in y
    'x',   #X vector
    'y',   #Y vector
    'valid'])


GlobalCalibration = namedtuple_with_defaults('GlobalCalibration', 
    ['umperpix', #Pixel size of the XTCAV camera
    'strstrength', #Strength parameter
    'rfampcalib', #Calibration of the RF amplitude
    'rfphasecalib', #Calibration of the RF phase
    'dumpe',        #Beam energy: dump config
    'dumpdisp',
    'valid'
    ], {'valid': 1})

      
ShotToShotParameters = namedtuple_with_defaults('ShotToShotParameters',
    ['ebeamcharge',  #ebeamcharge
    'dumpecharge',  #dumpecharge in C
    'xtcavrfamp',   #RF amplitude
    'xtcavrfphase', #RF phase
    'xrayenergy',         #Xrays energy in J
    'unixtime',
    'fiducial',
    'valid']
    )


ImageStatistics = namedtuple_with_defaults('ImageStatistics', 
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


PhysicalUnits = namedtuple_with_defaults('PhysicalUnits', 
    ['xfs',
    'yMeV',
    'xfsPerPix',
    'yMeVPerPix',
    'valid'])


Parameters = namedtuple_with_defaults('Parameters', 
    ['experiment', 'maxshots', 'run', 'validityrange', 
    'darkreferencepath', 'nb', 'groupsize', 'medianfilter', 
    'snrfilter', 'roiwaistthres', 'roiexpand', 'islandsplitmethod',
    'islandsplitpar1', 'islandsplitpar2', 'calpath', 'version'])

