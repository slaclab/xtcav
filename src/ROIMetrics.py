
import numpy as np

class ROIMetrics(object):
    

    def __init__(self, roiXN, roiX, roiYN, roiY):
        # Don't allow None arguments in class
        self.XN = roiXN  #Size of the image in X                           
        self.X = roiX    #Position of the first pixel in x
        self.YN = roiYN  #Size of the image in Y 
        self.Y = roiY    #Position of the first pixel in y
        
        self.valid = True
        if roiX is None:       
            warnings.warn_explicit('No XTCAV ROI info',UserWarning,'XTCAV',0)
            self.valid = False
            self.XN = 1024   #Size of the image in X                           
            self.X = 0         #Position of the first pixel in x
            self.YN = 1024   #Size of the image in Y 
            self.Y = 0         #Position of the first pixel in y

        self.x = roiX+np.arange(0, roiXN) #X vector
        self.y = roiY+np.arange(0, roiYN) #Y vector
        
        

    def Save(self,path):        
        constSave(self,path)

    @staticmethod    
    def Load(path):        
        return constLoad(path)
