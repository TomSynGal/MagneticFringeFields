###############################################################################
# Version 9.02 Params, Last Modified, 09/12/2021
###############################################################################
import numpy as np
###############################################################################

bjArray = np.array([0, 0.95, 0.90, 0.85], dtype = np.complex_)

lengthFactor = 1 # Alter the base unit of any function of length, 1 = 1m, 1E-03 = 1mm.

fieldFactor = 1 # Alter the base unit of any field, 1 = 1T, 1E-03 = 1mT.

zLim = 50.0 # The positive limit of Z in metres. (Lower limit is mirrored)

zLimPlot = 2.0 # Smaller limit for plots.

J = complex(0, 1)

###############################################################################

# Ue HFD5 instead of CSV. Structure based on magnet type and pole structure.

# Do Numpy gradient.