import numpy as np
import matplotlib.pyplot as plt
from V07MathematicaFunctions import bx, by, bzeta, ScalarPotential


# Let's pretend that everything works as intended.

# No questions asked! 

limit = 2.0 # Limit in the +/- z direction

iterations = 101

npts = 31

samples = 101

#bjArray = np.array([0, 0.95, 0.90, 0.85], dtype = np.complex_)
bjArray = np.array([0, 0.98, 0.95, 0.92], dtype = np.complex_) # BJ Array Changed in Newest Mathematica Notebook.

x = 0.5
y = 0.1
zRange = np.linspace(-limit,limit,iterations)

n = 2 # Multipole Number - Sextupole

bxData = np.zeros([iterations], dtype = np.complex_)
byData = np.zeros([iterations], dtype = np.complex_)
bzData = np.zeros([iterations], dtype = np.complex_)
scalarData = np.zeros([iterations], dtype = np.complex_)

# Create data loop
##############################################################################
for i in range(iterations):
    
    z = zRange[i]
    
    bxVal = bx(x, y, z, bjArray, n)
    bxData[i] = bxVal
    
    byVal = by (x, y, z, bjArray, n)
    byData[i]  = byVal
    
    bzVal = bzeta(x, y, z, bjArray, n)
    bzData[i] = bzVal
    
    scalarVal = ScalarPotential(x, y, z, bjArray, npts, samples, n)
    scalarData[i] = scalarVal
##############################################################################

scalarDiff = np.diff(scalarData, n=1)
scalarDiff = scalarDiff*18 # Some random scaling error (need to find the cause of this)
zRangeDiff = np.linspace(-limit,limit,(iterations-1))

BXfactor = 4
BYfactor = 9.7

plt.clf()
plt.title(r'Sextupole $\mathrm{\mathbb{Re}}$ Bx component as a function of z')
plt.plot(zRange, np.real(bxData), label = r'$\mathrm{\mathbb{Re}}$ Bx')
plt.xlabel('z (arb units)')
plt.ylabel('Output (arb units)')
plt.legend(loc = 2)
plt.show()

plt.clf()
plt.title(r'Sextupole $\mathrm{\mathbb{Re}}$ By component as a function of z')
plt.plot(zRange, np.real(byData), label = r'$\mathrm{\mathbb{Re}}$ By')
plt.xlabel('z (arb units)')
plt.ylabel('Output (arb units)')
plt.legend(loc = 2)
plt.show()

plt.clf()
plt.title(r'Sextupole $\mathrm{\mathbb{Re}}$ Bz component as a function of z')
plt.plot(zRange, np.real(bzData), label = r'$\mathrm{\mathbb{Re}}$ Bz')
plt.xlabel('z (arb units)')
plt.ylabel('Output (arb units)')
plt.legend(loc = 2)
plt.show()

plt.clf()
plt.title(r'Sextupole $\mathrm{\mathbb{Re}}$ scalar potential component as a function of z')
plt.plot(zRange, np.real(scalarData), label = r'$\mathrm{\mathbb{Re}}$ $\phi$')
plt.xlabel('z (arb units)')
plt.ylabel('Output (arb units)')
plt.legend(loc = 2)
plt.show()

plt.clf()
plt.title(r'Sextupole $\mathrm{\mathbb{Re}}$ $\partial_z$ scalar potential component as a function of z')
plt.plot(zRangeDiff, np.real(scalarDiff), label = r'$\partial_z$ ($\mathrm{\mathbb{Re}}$ $\phi$)')
plt.xlabel('z (arb units)')
plt.ylabel('Output (arb units)')
plt.legend(loc = 2)
plt.show()

plt.clf()
plt.title(r'Sextupole $\mathrm{\mathbb{Re}}$ Bx and $\partial_x$ ($\mathrm{\mathbb{Re}}$ $\phi$) component as a function of z')
plt.plot(zRange, np.real(bxData), label = r'$\mathrm{\mathbb{Re}}$ Bx')
plt.plot(zRange, (BXfactor*np.real(scalarData)), label = r'$\mathrm{\mathbb{Re}}$ $\phi$')
plt.xlabel('z (arb units)')
plt.ylabel('Output (arb units)')
plt.legend(loc = 2)
plt.show()

plt.clf()
plt.title(r'Sextupole $\mathrm{\mathbb{Re}}$ By and $\partial_y$ ($\mathrm{\mathbb{Re}}$ $\phi$) component as a function of z')
plt.plot(zRange, np.real(byData), label = r'$\mathrm{\mathbb{Re}}$ By')
plt.plot(zRange, (BYfactor*np.real(scalarData)), label = r'$\mathrm{\mathbb{Re}}$ $\phi$')
plt.xlabel('z (arb units)')
plt.ylabel('Output (arb units)')
plt.legend(loc = 2)
plt.show()

plt.clf()
plt.title(r'Sextupole $\mathrm{\mathbb{Re}}$ Bz and $\partial_z$ ($\mathrm{\mathbb{Re}}$ $\phi$) scalar potential component as a function of z')
plt.plot(zRange, np.real(bzData), label = r'$\mathrm{\mathbb{Re}}$ $\phi$')
plt.plot(zRangeDiff, np.real(scalarDiff), label = r'$\partial_z$ ($\mathrm{\mathbb{Re}}$ $\phi$)')
plt.xlabel('z (arb units)')
plt.ylabel('Output (arb units)')
plt.legend(loc = 2)
plt.show()
