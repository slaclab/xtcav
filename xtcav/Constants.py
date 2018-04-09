"""
Detector Names
"""
SRC = 'XrayTransportDiagnostic.0:Opal1000.0'
CALIB_GROUP = 'Xtcav::CalibV1'
ANALYSIS_VERSION = 'XTCAV_Analysis_Version'

EBEAM = 'EBeam'
GAS_DETECTOR = 'FEEGasDetEnergy'

ROI_SIZE_X = 'XTCAV_ROI_sizeX'
ROI_SIZE_Y = 'XTCAV_ROI_sizeY'
ROI_START_X = 'XTCAV_ROI_startX'
ROI_START_Y = 'XTCAV_ROI_startY'

"""
End Detector Names
"""

#Electron charge in coulombs
E_CHARGE=1.60217657e-19 

E_BEAM_CHARGE=5
XTCAV_RFAMP=20
XTCAV_RFPHASE=90
ENERGY_DETECTOR=0.2
DUMP_E_CHARGE=175E-12 #IN C

SNR_BORDER=100 #number of pixels near the border that can be considered to contain just noise
MIN_ROI_SIZE=3 #minimum number of pixels defining region of interest

DEFAULT_SPLIT_METHOD='scipyLabel'

DB_FILE_NAME = 'pedestals'
LOR_FILE_NAME = 'lasingoffreference'
