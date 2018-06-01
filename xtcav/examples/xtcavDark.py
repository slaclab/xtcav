from xtcav.DarkBackground import *
import numpy as np

dark_background = DarkBackground(
	experiment='amox23616', 
	run_number='104', 
	max_shots=500)

# Debugging
# print np.sum(dark_background.image)
