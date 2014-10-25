from xtcav.GenerateLasingOffReference import *

GLOC=GenerateLasingOffReference();
GLOC.experiment='amoc8114'
GLOC.runs='86'
GLOC.maxshots=401
GLOC.nb=1
GLOC.groupsize=5
GLOC.SetValidityRange(86,91)

GLOC.Generate();
