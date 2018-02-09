import numpy as np
import warnings

class ROIMetrics(object):
    

    def __init__(self, roiXN, roiX, roiYN, roiY, x=None, y=None):
        # Don't allow None arguments in class
        self.xN = roiXN  #Size of the image in X                           
        self.x0 = roiX    #Position of the first pixel in x
        self.yN = roiYN  #Size of the image in Y 
        self.y0 = roiY    #Position of the first pixel in y
        
        self.valid = True
        if roiX is None:       
            warnings.warn_explicit('No XTCAV ROI info',UserWarning,'XTCAV',0)
            self.valid = False
            self.xN = 1024   #Size of the image in X                           
            self.x0 = 0         #Position of the first pixel in x
            self.yN = 1024   #Size of the image in Y 
            self.y0 = 0         #Position of the first pixel in y

        self.x = x if x is not None else roiX+np.arange(0, roiXN) #X vector
        self.y = y if y is not None else roiY+np.arange(0, roiYN) #Y vector
        
        

    def Save(self,path):        
        constSave(vars(self),path)

    @staticmethod    
    def Load(path):        
        return constLoad(path)

class GlobalCalibration(object):
    

    def __init__(self, umperpix=None, strstrength=None, rfampcalib=None, rfphasecalib=None, dumpe=None, dumpdisp=None):
        self.umperpix=umperpix #Pixel size of the XTCAV camera
        self.strstrength=strstrength  #Strength parameter
        self.rfampcalib=rfampcalib    #Calibration of the RF amplitude
        self.rfphasecalib=rfphasecalib    #Calibration of the RF phase
        self.dumpe=dumpe                 #Beam energy: dump config
        self.dumpdisp=dumpdisp    

        self.valid = True
        for key, value in vars(self).iteritems():
            if value is None:
                warnings.warn_explicit('No XTCAV Calibration for epics variable' + key,UserWarning,'XTCAV',0)
                self.valid = False
        

    def Save(self,path):        
        constSave(vars(self),path)

    @staticmethod    
    def Load(path):        
        return constLoad(path)
