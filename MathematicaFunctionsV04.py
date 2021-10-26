import numpy as np
from scipy.special import factorial
from mpmath import polylog

#Make the complex number (0+i).
x1 = 0
y1 = 1
j = complex(x1, y1)


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


def bx(x, y, z, bjArray, n):
    
    output = (1/np.sqrt(2))*(bu(x, y, z, bjArray, n) + bv(x, y, z, bjArray, n))
    
    return output


def by(x, y, z, bjArray, n):
    
    output = (-j/np.sqrt(2))*(bu(x, y, z, bjArray, n) - bv(x, y, z, bjArray, n))
    
    return output
