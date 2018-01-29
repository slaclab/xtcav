
import numpy

class Parameters(object):
    def __init__(self):
        self.averagedProfiles=[]
        self.runs=numpy.array([],dtype=int)
        self.n=0
        self.parameters=0
        
    def Save(self,path):        
        constSave(self,path)

    @staticmethod    
    def Load(path):        
        return constLoad(path)
