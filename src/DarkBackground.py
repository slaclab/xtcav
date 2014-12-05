
import numpy
from xtcav.Constants import Load as constLoad
from xtcav.Constants import Save as constSave

class DarkBackground(object):
    def __init__(self):
        self.image=[]
        self.ROI=[]
        self.runs=numpy.array([],dtype=int)
        self.n=0
        
    def Save(self,path):        
        constSave(self,path)
        
    @staticmethod    
    def Load(path):        
        return constLoad(path)
