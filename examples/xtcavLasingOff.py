from xtcav.GenerateLasingOffReference import *

GLOC=GenerateLasingOffReference();
GLOC.experiment='amox23616'
GLOC.runs='131'
GLOC.maxshots=401
GLOC.nb=1
GLOC.groupsize=5
GLOC.SetValidityRange(86,91)

GLOC.Generate();
