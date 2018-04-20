from xtcav.DarkBackground import *
import numpy as np

dark_background = DarkBackground(
	experiment='diamcc14', 
	run_number='573', 
	maxshots=500,
	validityrange=(572,'end'))

# Debugging
# print np.sum(dark_background.image)
