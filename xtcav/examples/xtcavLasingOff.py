from xtcav.LasingOffReference import *

lasing_off_reference = LasingOffReference(
	experiment='amox23616',
	run_number='131',
	maxshots=100,
	num_bunches=1,
	groupsize=5)
lasing_off_reference.Generate()
print np.sum(lasing_off_reference.averaged_profiles.eCOMslice[0][0])
print "hi"
