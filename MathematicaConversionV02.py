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

x=1
y=1
z=1

#General Definitions

u = (1/np.sqrt(2))*complex(x,y)

v = (1/np.sqrt(2))*complex(x,-y)

zeta = np.sqrt(2)*z

def hj(bj1, u, v):
    
    output = (u/bj1)+(bj1*v)
    
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
print('The value of bj is...')
print(bj)


Product1 = (np.square(bj[[0]])-np.square(bj[[1]]))

cj = np.array([(np.power(-j, n)*np.power(bj[[0]], (n-1)))/(2*Product1),
               (np.power(-j, n)*np.power(bj[[1]], (n-1)))/-(2*Product1)])
print('The value of cj is...')
print(cj)



testt = zeta+ (j*(hj(1.05, u, v)))
print('here')
print(testt)

'''
def NewFunction(eta, n1):
    
    output = f(eta, n1)
'''

u = (1/np.sqrt(2))*complex(x,y)

print(u)
