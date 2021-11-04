##############################################################################
import numpy as np
import matplotlib.pyplot as plt
from V05MathematicaFunctions import G, bx, by, IntegrateTrapz
from V05MathematicaFunctions import NumericalFourierFit, FFT, FFTFit
##############################################################################


#Make the complex number (0+i).
x1 = 0
y1 = 1
j = complex(x1, y1)


# General Definitions.


data = np.zeros([5,200])

for i in range(5):
    for k in range (200):
        data[i][k] = G((k-100), i)

xdataRange = list(range(-100,100))
division = 2
xdataRange = [x/2 for x in xdataRange]

data = np.vstack((data, xdataRange))

test = data[1,:]
plt.title('General Definitions')
plt.plot(data[5,:], data[0,:], label = 'n=0')
plt.plot(data[5,:], data[1,:], label = 'n=1')
plt.plot(data[5,:], data[2,:], label = 'n=2')
plt.plot(data[5,:], data[3,:], label = 'n=3')
plt.plot(data[5,:], data[4,:], label = 'n=4')
plt.xlabel('z (arb units)')
plt.ylabel('Gradient Function')
plt.legend()
plt.show()

##############################################################################

# All Multipoles

bjArray = np.array([0, 0.95, 0.90, 0.85], dtype = np.complex_) # bj elements upto an octupole. # These elements are unitless.
# bjArray[[0]] element to be calculated by the bj function and is different for each multipole.

##############################################################################

# Quadrupole

n = 1 # Multipole number = 1 for the quadrupole. No units.

Zvalues = np.linspace(-3, 3, 6)
horizRange = np.linspace(-0.1, 0.1, 17)


# Graph 1
data1 = np.zeros([6,3,17])   
for i in range(6):
    for k in range (17):
        
        data1[i,0,:] = Zvalues[i]
        data1[i,1,k] = horizRange[k]
        mpc = bx(0, data1[i,1,k], data1[i,0,k], bjArray, 1)
        mpc = mpc.real
        data1[i,2,k] = mpc
        

plt.title('Quadrupole Bx Plot')
plt.plot(data1[0,1,:], data1[0,2,:], label = 'z=-3')
plt.plot(data1[1,1,:], data1[1,2,:], label = 'z=-1.8')
plt.plot(data1[2,1,:], data1[2,2,:], label = 'z=-0.6')
plt.plot(data1[3,1,:], data1[3,2,:], label = 'z=0.6')
plt.plot(data1[4,1,:], data1[4,2,:], label = 'z=1.8')
plt.plot(data1[5,1,:], data1[5,2,:], label = 'z=3')
plt.xlabel('y (arb units)')
plt.ylabel('Bx (arb units)')
plt.legend()
plt.show()

# Graph 2
data2 = np.zeros([6,3,17])   
for i in range(6):
    for k in range (17):
        
        data2[i,0,:] = Zvalues[i]
        data2[i,1,k] = horizRange[k]
        mpc = by(data2[i,1,k], 0, data2[i,0,k], bjArray, 1)
        #mpc = by(data2[i,1,k], 0, data2[i,0,k], 1)
        mpc = mpc.real
        data2[i,2,k] = mpc
        

plt.title('Quadrupole By Plot')
plt.plot(data2[0,1,:], data2[0,2,:], label = 'z=-3')
plt.plot(data2[1,1,:], data2[1,2,:], label = 'z=-1.8')
plt.plot(data2[2,1,:], data2[2,2,:], label = 'z=-0.6')
plt.plot(data2[3,1,:], data2[3,2,:], label = 'z=0.6')
plt.plot(data2[4,1,:], data2[4,2,:], label = 'z=1.8')
plt.plot(data2[5,1,:], data2[5,2,:], label = 'z=3')
plt.xlabel('x (arb units)')
plt.ylabel('By (arb units)')
plt.legend()
plt.show()

##############################################################################

# Sextupole

n = 2 # Multipole number = 2 for the sextupole. No units.

# Graph 3
data3 = np.zeros([6,3,17])   
for i in range(6):
    for k in range (17):
        
        data3[i,0,:] = Zvalues[i]
        data3[i,1,k] = horizRange[k]
        mpc = bx(0.1, data3[i,1,k], data3[i,0,k], bjArray, 2)
        mpc = mpc.real
        data3[i,2,k] = mpc
        

plt.title('Sextupole Bx Plot')
plt.plot(data3[0,1,:], data3[0,2,:], label = 'z=-3')
plt.plot(data3[1,1,:], data3[1,2,:], label = 'z=-1.8')
plt.plot(data3[2,1,:], data3[2,2,:], label = 'z=-0.6')
plt.plot(data3[3,1,:], data3[3,2,:], label = 'z=0.6')
plt.plot(data3[4,1,:], data3[4,2,:], label = 'z=1.8')
plt.plot(data3[5,1,:], data3[5,2,:], label = 'z=3')
plt.xlabel('y (arb units)')
plt.ylabel('Bx (arb units)')
plt.legend()
plt.show()

# Graph 4
data4 = np.zeros([6,3,17])   
for i in range(6):
    for k in range (17):
        
        data4[i,0,:] = Zvalues[i]
        data4[i,1,k] = horizRange[k]
        mpc = by(data4[i,1,k], 0, data4[i,0,k], bjArray, 2)
        #mpc = by(data2[i,1,k], 0, data2[i,0,k], 1)
        mpc = mpc.real
        data4[i,2,k] = mpc
        

plt.title('Sextupole By Plot')
plt.plot(data4[0,1,:], data4[0,2,:], label = 'z=-3')
plt.plot(data4[1,1,:], data4[1,2,:], label = 'z=-1.8')
plt.plot(data4[2,1,:], data4[2,2,:], label = 'z=-0.6')
plt.plot(data4[3,1,:], data4[3,2,:], label = 'z=0.6')
plt.plot(data4[4,1,:], data4[4,2,:], label = 'z=1.8')
plt.plot(data4[5,1,:], data4[5,2,:], label = 'z=3')
plt.xlabel('x (arb units)')
plt.ylabel('By (arb units)')
plt.legend()
plt.show()

##############################################################################

# Octupole

n = 3 # Multipole number = 2 for the sextupole. No units.

# Graph 5
data5 = np.zeros([6,3,17])   
for i in range(6):
    for k in range (17):
        
        data5[i,0,:] = Zvalues[i]
        data5[i,1,k] = horizRange[k]
        mpc = bx(0.1, data5[i,1,k], data5[i,0,k], bjArray, 3)
        mpc = mpc.real
        data5[i,2,k] = mpc
        

plt.title('Octupole Bx Plot')
plt.plot(data5[0,1,:], data5[0,2,:], label = 'z=-3')
plt.plot(data5[1,1,:], data5[1,2,:], label = 'z=-1.8')
plt.plot(data5[2,1,:], data5[2,2,:], label = 'z=-0.6')
plt.plot(data5[3,1,:], data5[3,2,:], label = 'z=0.6')
plt.plot(data5[4,1,:], data5[4,2,:], label = 'z=1.8')
plt.plot(data5[5,1,:], data5[5,2,:], label = 'z=3')
plt.xlabel('y (arb units)')
plt.ylabel('Bx (arb units)')
plt.legend()
plt.show()

# Graph 6
data6 = np.zeros([6,3,17])   
for i in range(6):
    for k in range (17):
        
        data6[i,0,:] = Zvalues[i]
        data6[i,1,k] = horizRange[k]
        mpc = by(data6[i,1,k], 0, data6[i,0,k], bjArray, 3)
        #mpc = by(data2[i,1,k], 0, data2[i,0,k], 1)
        mpc = mpc.real
        data6[i,2,k] = mpc
        

plt.title('Octupole By Plot')
plt.plot(data6[0,1,:], data6[0,2,:], label = 'z=-3')
plt.plot(data6[1,1,:], data6[1,2,:], label = 'z=-1.8')
plt.plot(data6[2,1,:], data6[2,2,:], label = 'z=-0.6')
plt.plot(data6[3,1,:], data6[3,2,:], label = 'z=0.6')
plt.plot(data6[4,1,:], data6[4,2,:], label = 'z=1.8')
plt.plot(data6[5,1,:], data6[5,2,:], label = 'z=3')
plt.xlabel('x (arb units)')
plt.ylabel('By (arb units)')
plt.legend()
plt.show()

##############################################################################

# Fit a Fourier series via numerical integration using Numpy Trapz.

# First calculate the Fourier coefficients.

hRange, data7 = IntegrateTrapz(21, 101, 2)
absoluteData7 = np.absolute(data7)
data7 = np.log(absoluteData7)

plt.clf()
plt.title('Fourier Coefficient Values for Numerical Integration')
fig = plt.plot(hRange, data7, 'ok')
plt.xlabel('Harmonic Number')
plt.ylabel('Ln(Abs( Fourier Coefficient))')
plt.show()

# Now call upon a function to fit the Fourier series.

# Define an imaginary part.
imagpt = 2.3j

# Generate standard Enge roll-off G for comparison.

data8 = np.zeros([100])
data11 = np.zeros([100])
for k in range (100):
    
    RealandImagData = (G((k-50+imagpt), 2))
    data8[k] = np.real(RealandImagData)
    data11[k] = np.imag(RealandImagData)
xdataRange8 = np.linspace(-50,50,100)

data8 = np.vstack((data8, xdataRange8))
data11 = np.vstack((data11, xdataRange8))

data9 = np.zeros([100])
data12 = np.zeros([100])
for i in range(100):
    datapoint= (NumericalFourierFit((i-50+imagpt), 21, 100, 2))
    data9[i] += np.real(datapoint)
    data12[i] += np.imag(datapoint)
xdataRange9 = np.linspace(-50,50,100)

data9 = np.vstack((data9, xdataRange9))
data12 = np.vstack((data12, xdataRange9))


# Plot the four graphs.

plt.title('Numerical Real Fourier Fit Comparison')
plt.plot(data8[1,:], data8[0,:], label = 'Output of G, n = 2')
plt.plot(data9[1,:], data9[0,:], label = 'Output of Fourier Fit, n = 2')
plt.xlabel('z (arb units)')
plt.ylabel('Re(G)')
plt.legend()
plt.show()

data10 = data9 - data8
plt.title('Numerical Real Fourier Fit Difference')
plt.plot(data9[1,:], data10[0,:])
plt.xlabel('z (arb units)')
plt.ylabel('Residual of Re(G) Fit')
plt.show()


plt.title('Numerical Imaginary Fourier Fit Comparison')
plt.plot(data11[1,:], data11[0,:], label = 'Output of G, n = 2')
plt.plot(data12[1,:], data12[0,:], label = 'Output of Fourier Fit, n = 2')
plt.xlabel('z (arb units)')
plt.ylabel('Im(G)')
plt.legend()
plt.show()

data13 = data12 - data11
plt.title('Numerical Imaginary Fourier Fit Difference')
plt.plot(data12[1,:], data13[0,:])
plt.xlabel('z (arb units)')
plt.ylabel('Residual of Im(G) Fit')
plt.show()

##############################################################################

# Fit a Fourier series via numerical integration using Numpy Trapz.

# First calculate the Fourier coefficients.

n = 2 # Multipole number.

npts = 31

fftData = FFT(npts, n)

fftDataGraph = fftData[1:(npts-1)]
fftDataGraph = np.absolute(fftDataGraph)
fftDataGraph = np.log(fftDataGraph)
fftDataGraph = fftDataGraph[fftDataGraph > -20]
#zRange = np.linspace(1, (2*(len(fftDataGraph))-1), len(fftDataGraph))
zRange = np.linspace(1, int(npts-2), int((npts-1)/2))

plt.clf()
plt.title('Fourier Coefficient Values for the Fast Fourier Transform')
fig = plt.plot(zRange, fftDataGraph, 'or')
plt.xlabel('Harmonic Number')
plt.ylabel('Ln(Abs( Fourier Coefficient))')
plt.show()

plt.clf()
plt.title('Comparison of Fourier Coefficients for Both Methods')
fig = plt.plot(hRange, data7, 'ok')
fig = plt.plot(zRange, fftDataGraph, 'or')
plt.xlabel('Harmonic Number')
plt.ylabel('Ln(Abs( Fourier Coefficient))')
plt.show()

# Now call upon a function to fit the FFT.

# Define an imaginary part.
imagpt = 2.3j

# Generate standard Enge roll-off G for comparison.

data13 = np.zeros([100])
data16 = np.zeros([100])
for k in range (100):
    
    RealandImagData = (G((k-50+imagpt), 2))
    data13[k] = np.real(RealandImagData)
    data16[k] = np.imag(RealandImagData)
xdataRange13 = np.linspace(-50,50,100)

data13 = np.vstack((data13, xdataRange13))
data16 = np.vstack((data16, xdataRange13))

data14 = np.zeros([100])
data17 = np.zeros([100])
for i in range(100):
    datapoint= (FFTFit((i-50+imagpt), 31, 100, 2))
    data14[i] += np.real(datapoint)
    data17[i] += np.imag(datapoint)
xdataRange14 = np.linspace(-50,50,100)

data14 = np.vstack((data14, xdataRange14))
data17 = np.vstack((data17, xdataRange14))


# Plot the four graphs.

plt.title('Fast Fourier Transform Real Fit Comparison')
plt.plot(data13[1,:], data13[0,:], label = 'Output of G, n = 2')
plt.plot(data14[1,:], data14[0,:], label = 'Output of Fourier Fit, n = 2')
plt.xlabel('z (arb units)')
plt.ylabel('Re(G)')
plt.legend()
plt.show()

data15 = data14 - data13
plt.title('Fast Fourier Transform Real Fit Difference')
plt.plot(data14[1,:], data15[0,:])
plt.xlabel('z (arb units)')
plt.ylabel('Residual of Re(G) Fit')
plt.show()


plt.title('Fast Fourier Transform Imaginary Fit Comparison')
plt.plot(data16[1,:], data16[0,:], label = 'Output of G, n = 2')
plt.plot(data17[1,:], data17[0,:], label = 'Output of Fourier Fit, n = 2')
plt.xlabel('z (arb units)')
plt.ylabel('Im(G)')
plt.legend()
plt.show()

data18 = data17 - data16
plt.title('Fast Fourier Transform Imaginary Fit Difference')
plt.plot(data17[1,:], data18[0,:])
plt.xlabel('z (arb units)')
plt.ylabel('Residual of Im(G) Fit')
plt.show()
