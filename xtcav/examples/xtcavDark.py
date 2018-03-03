from xtcav.DarkBackground import *
import numpy as np

dark_background = DarkBackground(
	experiment='amox23616', 
	run_number='104', 
	maxshots=150)

dark_background.generate()

#print np.sum(dark_background.image)