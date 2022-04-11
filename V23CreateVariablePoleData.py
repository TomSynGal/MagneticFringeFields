# Version 23.01 Create Variable Pole Data, Last Modified, 10/04/2022
###############################################################################
import h5py
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as interpolate
from V23Functions import ScalarPotential, ScalarPotentialSpline
###############################################################################

with h5py.File('C:/Users/thoma/Documents/#Python Code/PhD/Year 1/HDF5/VariableRoll.h5', 'r') as hdf:
    data1 = hdf.get('coeff')
    data2 = hdf.get('roll')
    coeff=np.array(data1)
    roll=np.array(data2)

#zvals = np.linspace(-2,2,201)
zvals = np.linspace(-50,50,2525)

plt.clf()
plt.plot(zvals,roll[0,:])
plt.plot(zvals,roll[1,:])
plt.plot(zvals,roll[2,:])
plt.plot(zvals,roll[3,:])
plt.plot(zvals,roll[4,:])
plt.show()

###############################################################################

nFourier = 45

XYplaneLim = 2
ZplaneLim = 50
Res = 11
zRes = 51

Scalar3D = np.zeros((len(roll),zRes,Res,Res))
XYvals = np.linspace(-XYplaneLim, XYplaneLim, Res)
Zvals = np.linspace(-ZplaneLim, ZplaneLim, zRes)

for h in range(len(roll)):
    
    print('Run, ',(h+1),' of ',len(roll))
    
    t, c, k = interpolate.splrep(zvals, roll[h,:], s=0, k=5)
    spline = interpolate.BSpline(t, c, k, extrapolate=False)
    
    for i in range(zRes):
        for j in range(Res):
            for k in range(Res):
                
                dataPoint = ScalarPotentialSpline(XYvals[j], XYvals[k], Zvals[i], nFourier, spline)
                Scalar3D[h,i,j,k] = dataPoint

ReshapeScalar3D = Scalar3D.reshape(len(roll),(zRes*Res*Res)) # Reshape to 2D

with h5py.File('C:/Users/thoma/Documents/#Python Code/PhD/Year 1/HDF5/VariablePole.h5', 'w') as hdf:
    hdf.create_dataset('coeff', data=coeff)
    hdf.create_dataset('roll', data=roll)
    hdf.create_dataset('scalar', data=ReshapeScalar3D)

###############################################################################

# Check that all is equal when the data is reformatted to 4D.

Scalar3D_1 = ReshapeScalar3D.reshape(len(roll),zRes,Res,Res)

test1 = Scalar3D[0,26,:,:]

test2 = Scalar3D_1[0,26,:,:]

###############################################################################