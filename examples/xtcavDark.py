from xtcav.GenerateDarkBackground import *

dark_background = generateDarkBackground(
	experiment='amox23616', 
	run_number='131', 
	maxshots=150, 
	validityrange=[85,109])
print dark_background.n
