from xtcav.DarkBackgroundReference import *
import numpy as np

dark_background = DarkBackgroundReference(
	experiment='cxin7316', 
	run_number='85', 
	max_shots=500)

# Debugging
# print np.sum(dark_background.image)