###############################################################################
# Version 10.01 HDF5 Test, Last Modified, 14/12/2021
###############################################################################
import numpy as np
import h5py
###############################################################################

matrix1 = np.random.random(size = (1000,1000))
matrix2 = np.random.random(size = (10000,100))

with h5py.File(r'C:\Users\thoma\Documents\#Python Code\PhD\Year 1\MagneticFringeFields Code\Week 10\hdf5test.h5', 'w') as hdf:
    hdf.create_dataset('dataset1', data = matrix1)
    hdf.create_dataset('dataset2', data = matrix2)
    

matrix3 = np.random.random(size = (1000,1000))
matrix4 = np.random.random(size = (1000,1000))
matrix5 = np.random.random(size = (1000,1000))
matrix6 = np.random.random(size = (1000,1000))

with h5py.File(r'C:\Users\thoma\Documents\#Python Code\PhD\Year 1\MagneticFringeFields Code\Week 10\hdf5test2.h5', 'w') as hdf:
    G1 = hdf.create_group('Group1')
    G1.create_dataset('dataset3', data = matrix3)
    G1.create_dataset('dataset4', data = matrix4)
    
    G2 = hdf.create_group('Group2\SubGroup1')
    G2.create_dataset('dataset5', data = matrix5)
    
    G3 = hdf.create_group('Group2\SubGroup2')
    G3.create_dataset('dataset6', data = matrix6)

dipo0 = np.random.random(size = (1000,1000))
dipo1 = np.random.random(size = (1000,1000))

quad0 = np.random.random(size = (1000,1000))
quad1 = np.random.random(size = (1000,1000))
quad2 = np.random.random(size = (1000,1000))
quad3 = np.random.random(size = (1000,1000))

sext0 = np.random.random(size = (1000,1000))
sext1 = np.random.random(size = (1000,1000))
sext2 = np.random.random(size = (1000,1000))
sext3 = np.random.random(size = (1000,1000))
sext4 = np.random.random(size = (1000,1000))
sext5 = np.random.random(size = (1000,1000))

octu0 = np.random.random(size = (1000,1000))
octu1 = np.random.random(size = (1000,1000))
octu2 = np.random.random(size = (1000,1000))
octu3 = np.random.random(size = (1000,1000))
octu4 = np.random.random(size = (1000,1000))
octu5 = np.random.random(size = (1000,1000))
octu6 = np.random.random(size = (1000,1000))
octu7 = np.random.random(size = (1000,1000))

with h5py.File(r'C:\Users\thoma\Documents\#Python Code\PhD\Year 1\MagneticFringeFields Code\Week 10\PoleFaces3D.h5', 'w') as hdf:
    G1 = hdf.create_group('Dipole')
    G1.create_dataset('Pole0', data = dipo0)
    G1.create_dataset('Pole1', data = dipo1)
    
    G2 = hdf.create_group('Quadrupole')
    G2.create_dataset('Pole0', data = quad0)
    G2.create_dataset('Pole1', data = quad1)
    G2.create_dataset('Pole2', data = quad2)
    G2.create_dataset('Pole3', data = quad3)
    
    G3 = hdf.create_group('Sextupole')
    G3.create_dataset('Pole0', data = sext0)
    G3.create_dataset('Pole1', data = sext1)
    G3.create_dataset('Pole2', data = sext2)
    G3.create_dataset('Pole3', data = sext3)
    G3.create_dataset('Pole4', data = sext4)
    G3.create_dataset('Pole5', data = sext5)

    G4 = hdf.create_group('Octupole')
    G4.create_dataset('Pole0', data = octu0)
    G4.create_dataset('Pole1', data = octu1)
    G4.create_dataset('Pole2', data = octu2)
    G4.create_dataset('Pole3', data = octu3)
    G4.create_dataset('Pole4', data = octu4)
    G4.create_dataset('Pole5', data = octu5)
    G4.create_dataset('Pole6', data = octu6)
    G4.create_dataset('Pole7', data = octu7)
