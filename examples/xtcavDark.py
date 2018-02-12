from xtcav.GenerateDarkBackground import *
import numpy as np

dark_background = DarkBackground(
	experiment='amox23616', 
	run_number='131', 
	maxshots=150, 
	validityrange=(85,109))
dark_background.Generate()
print np.sum(dark_background.image)