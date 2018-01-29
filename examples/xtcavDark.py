from xtcav.GenerateDarkBackground import *

GDB=GenerateDarkBackground();

GDB.experiment='amox23616'
GDB.runs='11'
GDB.maxshots=150
GDB.SetValidityRange(85,109)

GDB.Generate();
