import numpy as np
import math
from scipy.special import factorial, expn
from mpmath import polylog

#Make the complex number (0+i).
x1 = 0
y1 = 1
j = complex(x1, y1)

lim = 50.0

lim2 = 50.0

lim3 = 50.0


def u(x, y): # Units of u are metres, as x & y are in metres.
    
    u = (1/np.sqrt(2))*complex(x,y)
    
    return u

# Returns the complex number "u" in the form (x + iy).


def v(x, y): # Units of v are metres, as x & y are in metres.
    
    v = (1/np.sqrt(2))*complex(x,-y)
    
    return v

# Returns the complex number "v" in the form (x - iy).


def zeta(z): # Units of zeta are in metres, as z is in metres.
    
    zeta = np.sqrt(2)*z
    
    return zeta

# Returns zeta, the product of the z-coordinate and root 2.


def bj(bjArray, n):
    
    bj = bjArray[0:(n+1)] #Get only the needed elements of bjArray.
    
    bjToSum = bj[1:]
    
    numerator = -1*np.power(j, (n+1))
    
    denominator = np.prod(bjToSum)
    
    bjElementOne = numerator/denominator
    
    bj[0] = bjElementOne
    
    return bj


def cj(bjArray, n):
    
    zeroArray1 = np.array([0.0, 0.0, 0.0, 0.0], dtype = np.complex_)
        
    output = zeroArray1[0:(n+1)]
    
    for i in range (n+1):
        
        bjVals = bj(bjArray, n)
        
        halfNumerator = (1/2)*(np.power(-j, n)*np.power(bjVals[[i]],(n-1)))
        
        zeroArray2 = np.array([0.0, 0.0, 0.0, 0.0])
        
        seriesElements = zeroArray2[0:(n+1)]
        
        for k in range(n+1):
            
            if k == i:
                
                seriesElements[k] = 1
                
            else:
                
                seriesElements[k] = (np.square(bjVals[i])-np.square(bjVals[k]))
        
        summation = np.prod(seriesElements)
        
        output[i] = (halfNumerator/summation)

    return output


def hj(x, y, bjArray, bjIndex, n): # Units of hj are metres, as u & v are in metres and
                     # all b parameters seem to be unitless variables.
    
    bjVals = bj(bjArray, n)
    
    output = (u(x,y)/bjVals[bjIndex])+(bjVals[bjIndex]*v(x,y))
    
    return output

# Returns the sum of u & v, both modified by the parameter bj.


def G(xi, n1):
    
    if xi == 0:
        
        function = 0.5
        
    else:
        
        firstTerm = -(factorial(n1)/np.power(xi, n1))*polylog(n1, -np.exp(xi))
        
        secondTerm = 0
        
        for i in range (n1):
            
            k = i+1
            
            formula = (factorial(n1)/(factorial(n1-k)*np.power(xi, k)))*polylog(k, -1.0)
            
            secondTerm += formula
        
        function = firstTerm + secondTerm
    
    return function


def f(eta, n1):
    
    output = np.power(eta, n1)*G(eta, n1)
    
    return output


def g(eta, n1):
    
    output = np.power((-1), (n1+1))*np.power(eta, n1)*G(eta, n1)
    
    return output


def bu(x, y, z, bjArray, n):
    
    bjVals = bj(bjArray, n)
    
    cjVals = cj(bjArray, n)
    
    output = 0
    
    for i in range(n+1):
        
        firstFunc = cjVals[i]*bjVals[i]*f((zeta(z) + j*hj(x, y, bjArray, i, n)), n)
        
        secondFunc = -cjVals[i]*bjVals[i]*g((zeta(z) - j*hj(x, y, bjArray, i, n)), n)
        
        output += j*(firstFunc + secondFunc)
    
    return output


def bv(x, y, z, bjArray, n):
    
    bjVals = bj(bjArray, n)
    
    cjVals = cj(bjArray, n)
    
    output = 0
    
    for i in range(n+1):
        
        firstFunc = (cjVals[i]/bjVals[i])*f((zeta(z) + j*hj(x, y, bjArray, i, n)), n)
        
        secondFunc = -(cjVals[i]/bjVals[i])*g((zeta(z) - j*hj(x, y, bjArray, i, n)), n)
        
        output += j*(firstFunc + secondFunc)
    
    return output



def bzeta(x, y, z, bjArray, n):
    
    cjVals = cj(bjArray, n)
    
    output = 0
    
    for i in range(n+1):
        
        firstFunc = (cjVals[i])*f((zeta(z) + j*hj(x, y, bjArray, i, n)), n)
        
        secondFunc = (cjVals[i])*g((zeta(z) - j*hj(x, y, bjArray, i, n)), n)
        
        output += firstFunc + secondFunc
    
    return output


def bx(x, y, z, bjArray, n):
    
    output = (1/np.sqrt(2))*(bu(x, y, z, bjArray, n) + bv(x, y, z, bjArray, n))
    
    return output


def by(x, y, z, bjArray, n):
    
    output = (-j/np.sqrt(2))*(bu(x, y, z, bjArray, n) - bv(x, y, z, bjArray, n))
    
    return output


def IntegrateTrapz(hMax, samples, n):
    
    hRange = np.linspace(1, hMax, int((hMax+1)/2))
    data = np.zeros([hRange.size], dtype = np.complex_)
    
    for k in range(hRange.size):
        
            zSamples = samples
            zRange = np.linspace(-50.0, 50.0, zSamples)
            zMax = np.ndarray.max(zRange)
            zMin = np.ndarray.min(zRange)
            omega0 = np.pi/(zMax-zMin)
            
            zData = np.zeros([zSamples], dtype = np.complex_)
            
            for i in range (zSamples):
                
                nh = hRange[k]
            
                function = np.cos(nh*omega0*(zRange[i] + zMin))*(2*G((zRange[i] + 1E-06j), n) - 1)
                
                zData[i] = function
            
            integralVal = (omega0/np.pi)*np.trapz(zData, zRange)
                
            data[k] = integralVal
    
    return hRange, data


def NumericalFourierFit(xi, hMax, samples, n):
    
    zSamples = samples
    zRange = np.linspace(-lim, lim, zSamples)
    zMax = np.ndarray.max(zRange)
    zMin = np.ndarray.min(zRange)
    omega0 = np.pi/(zMax-zMin)
    
    hRange, data = IntegrateTrapz(hMax, samples, n)
    #hRange = int(hRange)
    #print(hRange)
    #print(hRange.size)
    length = int((hMax+1)/2)
    
    summation = 0.5
    
    for i in range(length):
        
        hValue = hRange[i]
        #print(hValue)
        
        function = data[i]*np.cos(hValue*omega0*(xi+zMin))
        
        summation += function
    
    return summation


def FFT(npts, n):
    
    zSamples = npts+1 
    zRange = np.linspace(-lim, lim, zSamples)
    
    data = np.zeros([zSamples], dtype = np.complex_)
    
    for i in range(zSamples):
        
        formula = (2*G(zRange[i], n) -1)
        
        data[i] = formula

    dataTake = data[1:(len(data)-1)]
    
    reverse = np.flip(dataTake)
    
    join = np.concatenate((data, reverse))

    fitCoeffs1 = (np.sqrt(len(join)))*(-1/(np.sqrt((2*npts)-1)))*np.fft.ifft(join)
    
    return fitCoeffs1


def FFTFit(xi, npts, samples, n):
    
    zSamples = samples
    zRange = np.linspace(-lim, lim, zSamples)
    zMax = np.ndarray.max(zRange)
    zMin = np.ndarray.min(zRange)
    omega0 = np.pi/(zMax-zMin)
    
    data = FFT(npts, n)
    
    length = int((npts-1)/2)
    
    nptsRange = np.linspace(1, int(npts-2), length)
    
    summation = 0.5
    
    for i in range(length):
        
        nptsValue = nptsRange[i]
        
        function = data[((2*i)+1)]*np.cos(nptsValue*omega0*(xi+zMin))
        
        summation += function
    
    return summation



def ScalarPotentialIntegral(xi, omega, samples, n):
    
    zSamples = samples
    zRange = np.linspace(-lim2, lim2, zSamples)
    zMin = np.ndarray.min(zRange)
    A = omega*zMin
    
    summation = 0
    
    for k in range (n+1):
        
        integral = (math.factorial(n)/math.factorial(n-k))*(np.power(xi, (n-k))/np.power(omega, (k+1)))*np.sin((omega*xi)+(0.5*k*np.pi)+A)
        
        summation += integral
    
    return summation



def fTilde(xi, npts, samples, n):
    
    zSamples = samples
    zRange = np.linspace(-lim3, lim3, zSamples)
    zMax = np.ndarray.max(zRange)
    zMin = np.ndarray.min(zRange)
    omega0 = np.pi/(zMax-zMin)
    
    coeffs = FFT(npts, n)
    
    length = int((npts-1)/2)
  
    summation = ((np.power((xi+zMin),(n+1)))/(2*(n+1)))
    
    for i in range(length):
        
        coeff = coeffs[((2*i)+1)]
        
        nh = (i+1)
        
        omega = ((2*nh)-1)*omega0
        
        #omega = omega0
        
        func = coeff*ScalarPotentialIntegral(xi, omega, samples, n)
        
        summation += func
    
    return summation



def ScalarPotential(x, y, z, bjArray, npts, samples, n):
    
    cjVals = cj(bjArray, n)
    
    output = 0
    
    for i in range(n+1):
        
        firstFunc = (cjVals[i])*fTilde((zeta(z) + j*hj(x, y, bjArray, i, n)), npts, samples, n)
        
        secondFunc = (cjVals[i])*(np.power(-1,(n+1)))*fTilde((zeta(z) - j*hj(x, y, bjArray, i, n)), npts, samples, n)
                
        output += firstFunc + secondFunc
    
    return output