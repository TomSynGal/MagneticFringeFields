###############################################################################
# Version 23.01 Params, Last Modified, 10/04/2022
###############################################################################
import numpy as np
###############################################################################

bjArray = np.array([0, 0.95, 0.90, 0.85], dtype = np.complex_)

LF = lengthFactor = 1 #E-03 # Alter the length base unit for all applications, 1 = 1m, 1E-03 = 1mm.

FF = fieldFactor = 1 #E-03 # Alter the field base unit for all applications, 1 = 1T, 1E-03 = 1mT.

zLim = 50.0 #15 #50.0 # The positive limit of Z in the base units stated above. (Lower limit is mirrored).

zLimVariable = 2.0 # Limit for the new variable created data.

# ZLim has always been 50, changed to 15 for the theoretical new ML data

zLimPlot = 50.0 # A smaller z limit for plotting in the base units stated above.

zLimVarG = 2.0

J = complex(0, 1) # Define the imaginary number j as J

###############################################################################

### PROGRESS REPORT ###

###############################################################################

