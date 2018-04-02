from xtcav.DarkBackground import *
import numpy as np

dark_background = DarkBackground(
	experiment='amox23616', 
	run_number='12', 
	maxshots=150,
	validityrange=(12,79))

# Debugging
# print np.sum(dark_background.image)
