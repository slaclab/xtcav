from xtcav.GenerateDarkBackground import *
import numpy as np

dark_background = generateDarkBackground(
	experiment='amoc8114', 
	run_number='110', 
	maxshots=150, 
	validityrange=[85,109])
print np.sum(dark_background.image)
