###############################################################################
# Version 9.01 Matlab3dData, Last Modified, 09/12/2021
###############################################################################
import numpy as np
import matplotlib.pyplot as plt
from V09Functions import ScalarPotential
###############################################################################

do_Dipole = False

do_Quadrupole = False

do_Sextupole = False

do_Octupole = True

###############################################################################
fourierRes = 101 #201
nFourierCoeffs = 45 #45
depth = 0.75
depthRes = 21
scalarPotentialLevel = 0.2
xLim = 1.5
yLim = 1.5
###############################################################################

xVals = np.linspace(-xLim, xLim, fourierRes)
yVals = np.linspace(-yLim, yLim, fourierRes)
zVals = np.linspace(-depth, depth, depthRes)
X , Y = np.meshgrid(xVals,yVals)

###############################################################################

if do_Dipole:
    
    datasetDipole = np.zeros([depthRes,fourierRes,fourierRes])
    
    diPole0Data = np.empty((0,3), float)
    diPole1Data = np.empty((0,3), float)

    for i in range(depthRes):
        for j in range(fourierRes):
            print('Layer ',(i+1),', Run ',(j+1),' of ',fourierRes)
            for k in range(fourierRes):
                
                z = zVals[i]
                x = X[0,j]
                y = Y[k,0]
                
                dataPoint = ScalarPotential(x, y, z, nFourierCoeffs, 0)
                datasetDipole[i,j,k] = dataPoint
    
        plt.clf()
        D = plt.contour(X, Y, np.absolute(datasetDipole[i,:,:]), levels = [scalarPotentialLevel])
        plt.show()
        pole0 = D.collections[0].get_paths()[0]
        pole1 = D.collections[0].get_paths()[1]
        
        pole0Coords = pole0.vertices
        pole0Coords = np.c_[pole0Coords,np.full((len(pole0Coords)), zVals[i]) ]
        diPole0Data = np.append(diPole0Data, pole0Coords, axis=0)
        pole0Coords = pole0.vertices
        
        pole1Coords = pole1.vertices
        pole1Coords = np.c_[pole1Coords,np.full((len(pole1Coords)), zVals[i]) ]
        diPole1Data = np.append(diPole1Data, pole1Coords, axis=0)
        pole1Coords = pole1.vertices
        
    np.savetxt("dipole0.csv", diPole0Data, delimiter=",")
    np.savetxt("dipole1.csv", diPole1Data, delimiter=",")

###############################################################################

if do_Quadrupole:
    
    datasetQuadrupole = np.zeros([depthRes,fourierRes,fourierRes])
    
    quadPole0Data = np.empty((0,3), float)
    quadPole1Data = np.empty((0,3), float)
    quadPole2Data = np.empty((0,3), float)
    quadPole3Data = np.empty((0,3), float)

    for i in range(depthRes):
        for j in range(fourierRes):
            print('Layer ',(i+1),', Run ',(j+1),' of ',fourierRes)
            for k in range(fourierRes):
                
                z = zVals[i]
                x = X[0,j]
                y = Y[k,0]
                
                dataPoint = ScalarPotential(x, y, z, nFourierCoeffs, 1)
                datasetQuadrupole[i,j,k] = dataPoint
    
        plt.clf()
        Q = plt.contour(X, Y, np.absolute(datasetQuadrupole[i,:,:]), levels = [scalarPotentialLevel])
        plt.show()
        pole0 = Q.collections[0].get_paths()[0]
        pole1 = Q.collections[0].get_paths()[1]
        pole2 = Q.collections[0].get_paths()[2]
        pole3 = Q.collections[0].get_paths()[3]
        
        pole0Coords = pole0.vertices
        pole0Coords = np.c_[pole0Coords,np.full((len(pole0Coords)), zVals[i]) ]
        quadPole0Data = np.append(quadPole0Data, pole0Coords, axis=0)
        pole0Coords = pole0.vertices
        
        pole1Coords = pole1.vertices
        pole1Coords = np.c_[pole1Coords,np.full((len(pole1Coords)), zVals[i]) ]
        quadPole1Data = np.append(quadPole1Data, pole1Coords, axis=0)
        pole1Coords = pole1.vertices
        
        pole2Coords = pole2.vertices
        pole2Coords = np.c_[pole2Coords,np.full((len(pole2Coords)), zVals[i]) ]
        quadPole2Data = np.append(quadPole2Data, pole2Coords, axis=0)
        pole2Coords = pole2.vertices
        
        pole3Coords = pole3.vertices
        pole3Coords = np.c_[pole3Coords,np.full((len(pole3Coords)), zVals[i]) ]
        quadPole3Data = np.append(quadPole3Data, pole3Coords, axis=0)
        pole3Coords = pole3.vertices
        
    np.savetxt("quad0.csv", quadPole0Data, delimiter=",")
    np.savetxt("quad1.csv", quadPole1Data, delimiter=",")
    np.savetxt("quad2.csv", quadPole2Data, delimiter=",")
    np.savetxt("quad3.csv", quadPole3Data, delimiter=",")

###############################################################################

if do_Sextupole:
    
    datasetSextupole = np.zeros([depthRes,fourierRes,fourierRes])
    
    sextuPole0Data = np.empty((0,3), float)
    sextuPole1Data = np.empty((0,3), float)
    sextuPole2Data = np.empty((0,3), float)
    sextuPole3Data = np.empty((0,3), float)
    sextuPole4Data = np.empty((0,3), float)
    sextuPole5Data = np.empty((0,3), float)

    for i in range(depthRes):
        for j in range(fourierRes):
            print('Layer ',(i+1),', Run ',(j+1),' of ',fourierRes)
            for k in range(fourierRes):
                
                z = zVals[i]
                x = X[0,j]
                y = Y[k,0]
                
                dataPoint = ScalarPotential(x, y, z, nFourierCoeffs, 2)
                datasetSextupole[i,j,k] = dataPoint
    
        plt.clf()
        S = plt.contour(X, Y, np.absolute(datasetSextupole[i,:,:]), levels = [scalarPotentialLevel])
        plt.show()
        pole0 = S.collections[0].get_paths()[0]
        pole1 = S.collections[0].get_paths()[1]
        pole2 = S.collections[0].get_paths()[2]
        pole3 = S.collections[0].get_paths()[3]
        pole4 = S.collections[0].get_paths()[4]
        pole5 = S.collections[0].get_paths()[5]
        
        pole0Coords = pole0.vertices
        pole0Coords = np.c_[pole0Coords,np.full((len(pole0Coords)), zVals[i]) ]
        sextuPole0Data = np.append(sextuPole0Data, pole0Coords, axis=0)
        pole0Coords = pole0.vertices
        
        pole1Coords = pole1.vertices
        pole1Coords = np.c_[pole1Coords,np.full((len(pole1Coords)), zVals[i]) ]
        sextuPole1Data = np.append(sextuPole1Data, pole1Coords, axis=0)
        pole1Coords = pole1.vertices
        
        pole2Coords = pole2.vertices
        pole2Coords = np.c_[pole2Coords,np.full((len(pole2Coords)), zVals[i]) ]
        sextuPole2Data = np.append(sextuPole2Data, pole2Coords, axis=0)
        pole2Coords = pole2.vertices
        
        pole3Coords = pole3.vertices
        pole3Coords = np.c_[pole3Coords,np.full((len(pole3Coords)), zVals[i]) ]
        sextuPole3Data = np.append(sextuPole3Data, pole3Coords, axis=0)
        pole3Coords = pole3.vertices
        
        pole4Coords = pole4.vertices
        pole4Coords = np.c_[pole4Coords,np.full((len(pole4Coords)), zVals[i]) ]
        sextuPole4Data = np.append(sextuPole4Data, pole4Coords, axis=0)
        pole4Coords = pole4.vertices
        
        pole5Coords = pole5.vertices
        pole5Coords = np.c_[pole5Coords,np.full((len(pole5Coords)), zVals[i]) ]
        sextuPole5Data = np.append(sextuPole5Data, pole5Coords, axis=0)
        pole5Coords = pole5.vertices
        
    np.savetxt("sextu0.csv", sextuPole0Data, delimiter=",")
    np.savetxt("sextu1.csv", sextuPole1Data, delimiter=",")
    np.savetxt("sextu2.csv", sextuPole2Data, delimiter=",")
    np.savetxt("sextu3.csv", sextuPole3Data, delimiter=",")
    np.savetxt("sextu4.csv", sextuPole4Data, delimiter=",")
    np.savetxt("sextu5.csv", sextuPole5Data, delimiter=",")

###############################################################################

if do_Octupole:
    
    datasetOctupole = np.zeros([depthRes,fourierRes,fourierRes])
    
    octuPole0Data = np.empty((0,3), float)
    octuPole1Data = np.empty((0,3), float)
    octuPole2Data = np.empty((0,3), float)
    octuPole3Data = np.empty((0,3), float)
    octuPole4Data = np.empty((0,3), float)
    octuPole5Data = np.empty((0,3), float)
    octuPole6Data = np.empty((0,3), float)
    octuPole7Data = np.empty((0,3), float)

    for i in range(depthRes):
        for j in range(fourierRes):
            print('Layer ',(i+1),', Run ',(j+1),' of ',fourierRes)
            for k in range(fourierRes):
                
                z = zVals[i]
                x = X[0,j]
                y = Y[k,0]
                
                dataPoint = ScalarPotential(x, y, z, nFourierCoeffs, 3)
                datasetOctupole[i,j,k] = dataPoint
    
        plt.clf()
        O = plt.contour(X, Y, np.absolute(datasetOctupole[i,:,:]), levels = [scalarPotentialLevel])
        plt.show()
        pole0 = O.collections[0].get_paths()[0]
        pole1 = O.collections[0].get_paths()[1]
        pole2 = O.collections[0].get_paths()[2]
        pole3 = O.collections[0].get_paths()[3]
        pole4 = O.collections[0].get_paths()[4]
        pole5 = O.collections[0].get_paths()[5]
        pole6 = O.collections[0].get_paths()[6]
        pole7 = O.collections[0].get_paths()[7]
        
        pole0Coords = pole0.vertices
        pole0Coords = np.c_[pole0Coords,np.full((len(pole0Coords)), zVals[i]) ]
        octuPole0Data = np.append(octuPole0Data, pole0Coords, axis=0)
        pole0Coords = pole0.vertices
        
        pole1Coords = pole1.vertices
        pole1Coords = np.c_[pole1Coords,np.full((len(pole1Coords)), zVals[i]) ]
        octuPole1Data = np.append(octuPole1Data, pole1Coords, axis=0)
        pole1Coords = pole1.vertices
        
        pole2Coords = pole2.vertices
        pole2Coords = np.c_[pole2Coords,np.full((len(pole2Coords)), zVals[i]) ]
        octuPole2Data = np.append(octuPole2Data, pole2Coords, axis=0)
        pole2Coords = pole2.vertices
        
        pole3Coords = pole3.vertices
        pole3Coords = np.c_[pole3Coords,np.full((len(pole3Coords)), zVals[i]) ]
        octuPole3Data = np.append(octuPole3Data, pole3Coords, axis=0)
        pole3Coords = pole3.vertices
        
        pole4Coords = pole4.vertices
        pole4Coords = np.c_[pole4Coords,np.full((len(pole4Coords)), zVals[i]) ]
        octuPole4Data = np.append(octuPole4Data, pole4Coords, axis=0)
        pole4Coords = pole4.vertices
        
        pole5Coords = pole5.vertices
        pole5Coords = np.c_[pole5Coords,np.full((len(pole5Coords)), zVals[i]) ]
        octuPole5Data = np.append(octuPole5Data, pole5Coords, axis=0)
        pole5Coords = pole5.vertices
        
        pole6Coords = pole6.vertices
        pole6Coords = np.c_[pole6Coords,np.full((len(pole6Coords)), zVals[i]) ]
        octuPole6Data = np.append(octuPole6Data, pole6Coords, axis=0)
        pole6Coords = pole6.vertices
        
        pole7Coords = pole7.vertices
        pole7Coords = np.c_[pole7Coords,np.full((len(pole7Coords)), zVals[i]) ]
        octuPole7Data = np.append(octuPole7Data, pole7Coords, axis=0)
        pole7Coords = pole7.vertices
        
    np.savetxt("octu0.csv", octuPole0Data, delimiter=",")
    np.savetxt("octu1.csv", octuPole1Data, delimiter=",")
    np.savetxt("octu2.csv", octuPole2Data, delimiter=",")
    np.savetxt("octu3.csv", octuPole3Data, delimiter=",")
    np.savetxt("octu4.csv", octuPole4Data, delimiter=",")
    np.savetxt("octu5.csv", octuPole5Data, delimiter=",")
    np.savetxt("octu6.csv", octuPole6Data, delimiter=",")
    np.savetxt("octu7.csv", octuPole7Data, delimiter=",")

###############################################################################