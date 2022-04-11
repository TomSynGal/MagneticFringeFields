###############################################################################
# Version 23.01 Create Variable Data, Last Modified, 10/04/2022
###############################################################################
import h5py
import numpy as np
from V23Functions import VarG
###############################################################################

Large = False

###############################################################################

if Large:
    samples = 10000

else:
    samples = 1000
        
    
res = 2525 # This will be 101 when the 50 range is reduced to 2
#zVals = np.linspace(-2,2,res)
zVals = np.linspace(-50,50,res)

data = np.zeros([samples,(res+3)])

for i in range(samples):
    
    print('Run ',(i+1),' of ',samples)
    
    a = 1.5*(np.random.uniform(low=-1.0, high=1.0))
    b = 3*(np.random.uniform(low=-1.0, high=1.0))
    c = 8*abs((np.random.uniform(low=-1.0, high=1.0)))
    output = np.transpose(VarG(a,b,c,res))
    output = np.column_stack((a,b,c,output))
    data[i] = output

coeff = data[:,0:3]

roll = data[:,3:]

if Large:
    with h5py.File('C:/Users/thoma/Documents/#Python Code/PhD/Year 1/HDF5/VariableRollLarge.h5', 'w') as hdf:
        hdf.create_dataset('coeff', data=coeff)
        hdf.create_dataset('roll', data=roll)

else:
    with h5py.File('C:/Users/thoma/Documents/#Python Code/PhD/Year 1/HDF5/VariableRoll.h5', 'w') as hdf:
        hdf.create_dataset('coeff', data=coeff)
        hdf.create_dataset('roll', data=roll)
    