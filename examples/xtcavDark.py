from xtcav.GenerateDarkBackground import *

GDB=GenerateDarkBackground();

GDB.experiment='amoc8114'
GDB.runs='85'
GDB.maxshots=150
GDB.SetValidityRange(85,109)

GDB.Generate();
