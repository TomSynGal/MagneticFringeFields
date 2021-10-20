#Mathematica notebook to Python conversion.

##############################################################################
import math
import cmath
import numpy as np
import mpmath as mp
import scipy
import sympy as sp
import matplotlib.pyplot as plt
from scipy.special import factorial
from mpmath import polylog
##############################################################################

#Parameters

#x=1
#y=1
#z=1

#General Definitions

def u(x, y):
    
    u = (1/np.sqrt(2))*complex(x,y)
    
    return u


def v(x, y):
    
    v = (1/np.sqrt(2))*complex(x,-y)
    
    return v


def zeta(z):
    
    zeta = np.sqrt(2)*z
    
    return zeta


def hj(bj1, x, y):
    
    output = (u(x,y)/bj1)+(bj1*v(x,y))
    
    return output


def G(xi, n1):
    
    firstTerm = -(factorial(n1)/np.power(xi, n1))*polylog(n1, -np.exp(xi))
    
    secondTerm = 0
    
    for i in range (n1):
        
        j = i+1
        
        formula = (factorial(n1)/(factorial(n1-j)*np.power(xi, j)))*polylog(j, -1.0)
        
        secondTerm += formula
    
    function = firstTerm + secondTerm
    
    return function


def f(eta, n1):
    
    output = np.power(eta, n1)*G(eta, n1)
    
    return output


def g(eta, n1):
    
    output = np.power(-1, (n1+1))*np.power(eta, n1)*G(eta, n1)
    
    return output

#General Definitions Plot
##############################################################################

#Data creation

data = np.zeros([5,200])


for i in range(5):
    for j in range (200):
        
        data[i][j] = G((j-100), i)

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

#Quadrupole


#Make the complex number (0+i)
x1 = 0
y1 = 1
j = complex(x1, y1)

n = 1
b1 = 1
bj = np.array([b1, 0.95])
bj[[0]] = -(np.power(j, (n+1))/bj[[1]])

Product1 = (np.square(bj[[0]])-np.square(bj[[1]]))
cj = np.array([(np.power(-j, n)*np.power(bj[[0]], (n-1)))/(2*Product1),
               (np.power(-j, n)*np.power(bj[[1]], (n-1)))/-(2*Product1)])

print('value cj')
print(cj)
'''

def cj(n):
    
    for i in range (n+1):
        
        numerator = (np.power(-j, n)*np.power(bj[[i]],(n-1)))
        
        denominator = 1
        
        for k in range (n+1):
            
            if k==i:
                
                denominator = denominator
            
            else:
                
                denominator = denominator*2*(np.square(bj[[i]])-np.square(bj[[k]]))
        
        cj = np.array([n+1])
        
        cj[i] = (numerator/denominator)
        
    return cj
    
        
    

'''
def bu(x, y, z, n):
    
    output = 0
    
    for i in range(n+1):
        
        firstFunc = cj[i]*bj[i]*f((zeta(z) + j*hj(bj[i], u(x,y), v(x,y))), n)
        
        secondFunc = -cj[i]*bj[i]*g((zeta(z) - j*hj(bj[i], u(x,y), v(x,y))), n)
        
        output += j*(firstFunc + secondFunc)
    
    return output


def bv(x, y, z, n):
    
    output = 0
    
    for i in range(n+1):
        
        firstFunc = (cj[i]/bj[i])*f((zeta(z) + j*hj(bj[i], u(x,y), v(x,y))), n)
        
        secondFunc = -(cj[i]/bj[i])*g((zeta(z) - j*hj(bj[i], u(x,y), v(x,y))), n)
        
        output += j*(firstFunc + secondFunc)
    
    return output


def bx(x, y, z, n):
    
    output = (1/np.sqrt(2))*(bu(x,y,z,n) + bv(x,y,z,n))
    
    return output


def by(x, y, z, n):
    
    output = (-j/np.sqrt(2))*(bu(x,y,z,n) - bv(x,y,z,n))
    
    return output



Zvalues = np.linspace(-3, 3, 6)
Yvalues = np.linspace(-0.1, 0.1, 9)
data1 = np.zeros([6,3,9])
    
for i in range(6):
    for k in range (9):
        
        data1[i,0,:] = Zvalues[i]
        data1[i,1,k] = Yvalues[k]
        mpc = -bx(0, data1[i,1,k], data1[i,0,k], 1)
        mpc = mpc[0].imag
        data1[i,2,k] = mpc
        

plt.title('Bx Plot')
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

data2 = np.zeros([6,3,9])
    
for i in range(6):
    for k in range (9):
        
        data2[i,0,:] = Zvalues[i]
        data2[i,1,k] = Yvalues[k]
        mpc = by(0, data2[i,1,k], data2[i,0,k], 1)
        mpc = mpc[0].imag
        data2[i,2,k] = mpc
        

plt.title('By Plot')
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

'''
##############################################################################

#Sextupole

x1 = 0
y1 = 1
j = complex(x1, y1)

n = 2
b1 = 1
bj = np.array([b1, 0.95, 0.90])
number = (-np.power(j, (n+1)))*(1/(bj[[1]]*bj[[2]]))
bj[[[0]]]=(number.imag)
print(number)

print('numberrr')
print(number.imag)

bj1 =(1/(bj[[1]]*bj[[2]]))*(-np.power(j, (n+1)))

Product1 = (np.square(bj[[0]])-np.square(bj[[1]]))
cj = np.array([(np.power(-j, n)*np.power(bj[[0]], (n-1)))/(2*Product1),
               (np.power(-j, n)*np.power(bj[[1]], (n-1)))/-(2*Product1)])

##############################################################################

#Octupole

x1 = 0
y1 = 1
j = complex(x1, y1)
'''