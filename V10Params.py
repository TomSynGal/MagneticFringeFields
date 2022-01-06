###############################################################################
# Version 10.01 Params, Last Modified, 14/12/2021
###############################################################################
import numpy as np
###############################################################################

bjArray = np.array([0, 0.95, 0.90, 0.85], dtype = np.complex_)

LF = lengthFactor = 1E-03 # Alter the base unit of any function of length, 1 = 1m, 1E-03 = 1mm.

FF = fieldFactor = 1E-03 # Alter the base unit of any field, 1 = 1T, 1E-03 = 1mT.

zLim = 50.0 # The positive limit of Z in metres. (Lower limit is mirrored)

zLimPlot = 2.0 # Smaller limit for plots.

J = complex(0, 1)

###############################################################################

### PROGRESS REPORT ###

# Changed numpy.diff to numpy.gradient, changes z scale and makes the code less bulky but still a need for a scaling factor.

# Changed saved data to HDF5 format.

# Implemented length and field factors hopefully correctly now. #upto 2d scalar potential plots.

#Okay fixed scaling factor, now need to find the error that overstates the differential of the scalar potential in z.