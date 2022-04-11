###############################################################################
# Version 23.01 Functions, Last Modified, 10/04/2022
###############################################################################
import math
import numpy as np
from scipy.special import factorial
from mpmath import polylog
from V23Params import bjArray, J, zLim, zLimVarG, zLimVariable
###############################################################################


def u(x, y):
    
    u = (1/np.sqrt(2))*complex(x,y)
    
    return u
# Returns the complex number "u" in the form (x + iy).


def v(x, y):
    
    v = (1/np.sqrt(2))*complex(x,-y)
    
    return v
# Returns the complex number "v" in the form (x - iy).


def zeta(z):
    
    zeta = np.sqrt(2)*z
    
    return zeta


def bj(n):
    
    bj = bjArray[0:(n+1)]
    
    bjToSum = bj[1:]
    
    numerator = -1*np.power(J, (n+1))
    
    denominator = np.prod(bjToSum)
    
    bjElementOne = numerator/denominator
    
    bj[0] = bjElementOne
    
    return bj
# Returns all of bj as an array of length (n+1).


def cj(n):
    
    zeroArray1 = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0], dtype = np.complex_)
        
    output = zeroArray1[0:(n+1)]
    
    for i in range (n+1):
        
        bjVals = bj(n)
        
        halfNumerator = (1/2)*(np.power(-J, n)*np.power(bjVals[[i]],(n-1)))
        
        zeroArray2 = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        
        seriesElements = zeroArray2[0:(n+1)]
        
        for k in range(n+1):
            
            if k == i:
                
                seriesElements[k] = 1
                
            else:
                
                seriesElements[k] = (np.square(bjVals[i])-np.square(bjVals[k]))
        
        summation = np.prod(seriesElements)
        
        output[i] = (halfNumerator/summation)

    return output
# Returns all of cj as an array of length (n+1).


def hj(x, y, bjIndex, n):
    
    bjVals = bj(n)
    
    output = (u(x,y)/bjVals[bjIndex])+(bjVals[bjIndex]*v(x,y))
    
    return output
# Returns a numerical value of hj for a given index of the bj array. So there are differnet values of hj for different bj inputs.
# Returns the sum of u & v, both modified by the parameter bj.


def G(xi, n):
    
    xi = -xi # REVERSAL FEATURE
    
    if xi == 0:
        
        function = 0.5
        
    else:
        
        firstTerm = -(factorial(n)/np.power(xi, n))*polylog(n, -np.exp(xi))
        
        secondTerm = 0
        
        for i in range (n):
            
            k = i+1
            
            formula = (factorial(n)/(factorial(n-k)*np.power(xi, k)))*polylog(k, -1.0)
            
            secondTerm += formula
        
        function = firstTerm + secondTerm
    
    return function

def VarG(a, b, c, res):
    
    zVals = np.linspace(-zLimVarG, zLimVarG, res)
    output = np.zeros([res,1])
    
    for i in range(res):
        
        output[i,0] = 1/(1+np.exp((a*zVals[i])+(b*np.power(zVals[i],3))+(c*np.power(zVals[i],5))))
    
    return output


def f(eta, n):
    
    output = np.power(eta, n)*G(eta, n)
    
    return output


def g(eta, n):
    
    output = np.power(-1, (n+1))*np.power(eta, n)*G(eta, n)
    
    return output


def fSpline(eta, Spline):
    
    n = 1 # Quadrupole
    
    output = np.power(eta, n)*Spline(eta)
    
    return output


def gSpline(eta, Spline):
    
    n = 1 # Quadrupole
    
    output = np.power(-1, (n+1))*np.power(eta, n)*Spline(eta)
    
    return output


def bu(x, y, z, n):
    
    bjVals = bj(n)
    
    cjVals = cj(n)
    
    output = 0
    
    for i in range(n+1):
        
        firstFunc = cjVals[i]*bjVals[i]*f((zeta(z) + J*hj(x, y, i, n)), n)
        
        secondFunc = -cjVals[i]*bjVals[i]*g((zeta(z) - J*hj(x, y, i, n)), n)
        
        output += J*(firstFunc + secondFunc)
    
    return output


def bv(x, y, z, n):
    
    bjVals = bj(n)
    
    cjVals = cj(n)
    
    output = 0
    
    for i in range(n+1):
        
        firstFunc = (cjVals[i]/bjVals[i])*f((zeta(z) + J*hj(x, y, i, n)), n)
        
        secondFunc = -(cjVals[i]/bjVals[i])*g((zeta(z) - J*hj(x, y, i, n)), n)
        
        output += J*(firstFunc + secondFunc)
    
    return output


def bzeta(x, y, z, n):
    
    cjVals = cj(n)
    
    output = 0
    
    for i in range(n+1):
        
        firstFunc = (cjVals[i])*f((zeta(z) + J*hj(x, y, i, n)), n)
        
        secondFunc = (cjVals[i])*g((zeta(z) - J*hj(x, y, i, n)), n)
        
        output += firstFunc + secondFunc
    
    return output


def buSpline(x, y, z, Spline):
    
    n = 1 # Quadrupole
    
    bjVals = bj(n)
    
    cjVals = cj(n)
    
    output = 0
    
    for i in range(n+1):
        
        firstFunc = cjVals[i]*bjVals[i]*fSpline((zeta(z) + J*hj(x, y, i, n)), Spline)
        
        secondFunc = -cjVals[i]*bjVals[i]*gSpline((zeta(z) - J*hj(x, y, i, n)), Spline)
        
        output += J*(firstFunc + secondFunc)
    
    return output


def bvSpline(x, y, z, Spline):
    
    n = 1 # Quadrupole
    
    bjVals = bj(n)
    
    cjVals = cj(n)
    
    output = 0
    
    for i in range(n+1):
        
        firstFunc = (cjVals[i]/bjVals[i])*fSpline((zeta(z) + J*hj(x, y, i, n)), Spline)
        
        secondFunc = -(cjVals[i]/bjVals[i])*gSpline((zeta(z) - J*hj(x, y, i, n)), Spline)
        
        output += J*(firstFunc + secondFunc)
    
    return output


def bzetaSpline(x, y, z, Spline):
    
    n = 1 # Quadrupole
    
    cjVals = cj(n)
    
    output = 0
    
    for i in range(n+1):
        
        firstFunc = (cjVals[i])*fSpline((zeta(z) + J*hj(x, y, i, n)), Spline)
        
        secondFunc = (cjVals[i])*gSpline((zeta(z) - J*hj(x, y, i, n)), Spline)
        
        output += firstFunc + secondFunc
    
    return output


def bx(x, y, z, n):
    
    output = (1/np.sqrt(2))*(bu(x, y, z, n) + bv(x, y, z, n))
    
    return output


def by(x, y, z, n):
    
    output = (-J/np.sqrt(2))*(bu(x, y, z, n) - bv(x, y, z, n))
    
    return output


def bxSpline(x, y, z, Spline):
    
    output = (1/np.sqrt(2))*(buSpline(x, y, z, Spline) + bvSpline(x, y, z, Spline))
    
    return output


def bySpline(x, y, z, Spline):
    
    output = (-J/np.sqrt(2))*(buSpline(x, y, z, Spline) - bvSpline(x, y, z, Spline))
    
    return output


def IntegrateTrapz(nFourier, res, n):
    
    coeffNumbers = np.linspace(1, nFourier, int((nFourier+1)/2))
    coeffs = np.zeros([coeffNumbers.size], dtype = np.complex_)
    zRange = np.linspace(-zLim, zLim, res)
    zMax = zLim
    zMin = -zLim
    omega0 = np.pi/(zMax-zMin)
    
    for k in range(coeffNumbers.size):
            
            zData = np.zeros([res], dtype = np.complex_)
            
            for i in range (res):
            
                function = np.cos(coeffNumbers[k]*omega0*(zRange[i] + zMin))*(2*G((zRange[i] + 1E-06j), n) - 1)
                
                zData[i] = function
            
            integralVal = (omega0/np.pi)*np.trapz(zData, zRange)
                
            coeffs[k] = integralVal
    
    return coeffNumbers, coeffs


def NumericalFourierFit(xi, nFourier, samples, n):
    
    zRange = np.linspace(-zLim, zLim, samples)
    zMax = zLim
    zMin = -zLim
    omega0 = np.pi/(zMax-zMin)
    
    coeffNumbers, coeffs = IntegrateTrapz(nFourier, samples, n)
    
    length = int((nFourier+1)/2)
    summation = 0.5
    
    for i in range(length):
        
        function = coeffs[i]*np.cos(coeffNumbers[i]*omega0*(xi+zMin))
        
        summation += function
    
    return summation, zRange


def FFT(npts, n):
    
    zSamples = npts+1 
    zRange = np.linspace(-zLim, zLim, zSamples)
    
    data = np.zeros([zSamples], dtype = np.complex_)
    
    for i in range(zSamples):
        
        formula = (2*G(zRange[i], n) -1)
        
        data[i] = formula

    dataTake = data[1:(len(data)-1)]
    
    reverse = np.flip(dataTake)
    
    join = np.concatenate((data, reverse))

    fitCoeffs1 = (np.sqrt(len(join)))*(-1/(np.sqrt((2*npts)-1)))*np.fft.ifft(join)
    
    return fitCoeffs1


def FFTSpline(npts, Spline):
    
    zSamples = npts+1 
    zRange = np.linspace(-zLim, zLim, zSamples)
    #zRange = np.linspace(-zLimVariable, zLimVariable, zSamples)
    
    data = np.zeros([zSamples], dtype = np.complex_)
    
    for i in range(zSamples):
        
        formula = (2*Spline(zRange[i]) -1)
        
        data[i] = formula

    dataTake = data[1:(len(data)-1)]
    
    reverse = np.flip(dataTake)
    
    join = np.concatenate((data, reverse))

    fitCoeffs1 = (np.sqrt(len(join)))*(-1/(np.sqrt((2*npts)-1)))*np.fft.ifft(join)
    
    return fitCoeffs1


def FFTFit(xi, npts, n):
    
    zMax = zLim
    zMin = -zLim
    omega0 = np.pi/(zMax-zMin)
    
    data = FFT(npts, n)
    
    length = int((npts-1)/2)
    
    nptsRange = np.linspace(1, int(npts-2), length)

    summation = 0.5
    
    for i in range(length):
        
        function = data[((2*i)+1)]*np.cos(nptsRange[i]*omega0*(xi+zMin))
        
        summation += function
    
    return summation


def FFTFitSpline(xi, npts, Spline):
    
    zMax = zLim
    zMin = -zLim
    #zMax = zLimVariable
    #zMin = -zLimVariable
    omega0 = np.pi/(zMax-zMin)
    
    data = FFTSpline(npts, Spline)
    
    length = int((npts-1)/2)
    
    nptsRange = np.linspace(1, int(npts-2), length)

    summation = 0.5
    
    for i in range(length):
        
        function = data[((2*i)+1)]*np.cos(nptsRange[i]*omega0*(xi+zMin))
        
        summation += function
    
    return summation


def ScalarPotentialIntegral(xi, omega, n):

    zMin = -zLim
    A = omega*zMin
    
    summation = 0
    
    for k in range (n+1):
        
        integral = (math.factorial(n)/math.factorial(n-k))*(np.power(xi, (n-k))/np.power(omega, (k+1)))*np.sin((omega*xi)+(0.5*k*np.pi)+A)
        
        summation += integral
    
    return summation



def fTilde(xi, npts, n):

    zMax = zLim
    zMin = -zLim
    omega0 = np.pi/(zMax-zMin)
    
    coeffs = FFT(npts, n)
    
    length = int((npts-1)/2)
  
    summation = 0.5+((np.power((xi+zMin),(n+1)))/(2*(n+1)))
    
    for i in range(length):
        
        coeff = coeffs[((2*i)+1)]
        
        nh = (i+1)
        
        omega = ((2*nh)-1)*omega0
        
        func = coeff*ScalarPotentialIntegral(xi, omega, n)
        
        summation += func
    
    return summation


def fTildeSpline(xi, npts, Spline):
    
    n = 1 ### QUADRUPOLE
    zMax = zLim
    zMin = -zLim
    #zMax = zLimVariable
    #zMin = -zLimVariable
    omega0 = np.pi/(zMax-zMin)
    
    coeffs = FFTSpline(npts, Spline)
    
    length = int((npts-1)/2)
  
    summation = 0.5+((np.power((xi+zMin),(n+1)))/(2*(n+1)))
    
    for i in range(length):
        
        coeff = coeffs[((2*i)+1)]
        
        nh = (i+1)
        
        omega = ((2*nh)-1)*omega0
        
        func = coeff*ScalarPotentialIntegral(xi, omega, n)
        
        summation += func
    
    return summation



def ScalarPotential(x, y, z, npts, n):
    
    cjVals = cj(n)
    
    output = 0
    
    for i in range(n+1):
        
        firstFunc = (cjVals[i])*fTilde((zeta(z) + J*hj(x, y, i, n)), npts, n)
        
        secondFunc = (cjVals[i])*(np.power(-1,(n+1)))*fTilde((zeta(z) - J*hj(x, y, i, n)), npts, n)
                
        output += firstFunc + secondFunc
    
    return output


def ScalarPotentialSpline(x, y, z, npts, Spline):
    
    n = 1 ### QUADRUPOLE
    cjVals = cj(n)
    
    output = 0
    
    for i in range(n+1):
        
        firstFunc = (cjVals[i])*fTildeSpline((zeta(z) + J*hj(x, y, i, n)), npts, Spline)
        
        secondFunc = (cjVals[i])*(np.power(-1,(n+1)))*fTildeSpline((zeta(z) - J*hj(x, y, i, n)), npts, Spline)
                
        output += firstFunc + secondFunc
    
    return output

###############################################################################