
import numpy
from Constants import Load as constLoad
from Constants import Save as constSave
import copy
from Metrics import *

class DarkBackground(object):
    def __init__(self):
        self.image=[]
        self.ROI=None
        self.run=''
        self.n=0
        
    def Save(self,path): 
        # super hacky... allows us to save 
        instance = copy.deepcopy(self)
        if instance.ROI:
            instance.ROI = vars(instance.ROI)
            instance.ROI.pop('valid')
        constSave(instance,path)
        
    @staticmethod    
    def Load(path):        
        db = constLoad(path)
        db.ROI = ROIMetrics(db.ROI['xN'], db.ROI['x0'], db.ROI['yN'], db.ROI['y0'], x=db.ROI['x'], y=db.ROI['y'])
        return db
