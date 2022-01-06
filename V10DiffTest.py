###############################################################################
# Version 10.01 Diff Test, Last Modified, 14/12/2021
###############################################################################
import numpy as np
import matplotlib.pyplot as plt
from V10Functions import bzeta, ScalarPotential
from V10Params import zLimPlot
###############################################################################
    
fourierRes = 1001
nFourierCoeffs = 45
diffFactor = 177
gradFactor = 177

x = 0.5
y = 0.1
zRange = np.linspace(-zLimPlot, zLimPlot, fourierRes)
zRangeDiff = np.linspace((-zLimPlot+(1/(2*fourierRes))),(zLimPlot-(1/(2*fourierRes))),(fourierRes-1)) # generate the new z scale.

n = 2 # Multipole Number - Sextupole

bzData = np.zeros([fourierRes], dtype = np.complex_)

bzDataDiffScale = np.zeros([fourierRes-1], dtype = np.complex_)
scalarPotData = np.zeros([fourierRes], dtype = np.complex_)

for i in range(fourierRes):
     
    bzDataPoint = bzeta(x, y, zRange[i], n)
    bzData[i] = bzDataPoint
   
for i in range(fourierRes-1):
    
    bzDataPointDiff = bzeta(x, y, zRangeDiff[i], n)
    bzDataDiffScale[i] = bzDataPointDiff

scalarPotData = ScalarPotential(x, y, zRange, nFourierCoeffs, n) # Compute first for FourierRes
scalarDiff = np.diff(scalarPotData, n=1)
scalarDiff = scalarDiff*diffFactor # Differentiate numerically  which alters the z scale.

scalarGrad = np.gradient(scalarPotData)
scalarGrad = scalarGrad*gradFactor

plt.clf()
plt.title(r'Sextupole $\mathrm{\mathbb{Re}}$ Bz and $\partial_z$ ($\mathrm{\mathbb{Re}}$ $\phi$) scalar potential component as a function of z')
plt.plot(zRangeDiff, np.real(bzDataDiffScale), label = r'$\mathrm{\mathbb{Re}}$ Bz')
plt.plot(zRangeDiff, np.real(scalarDiff), label = r'$\partial_z$ ($\mathrm{\mathbb{Re}}$ $\phi$)')
plt.xlabel(r'$z$ ($mm$)')
plt.ylabel(r'$\Re$ $B_{z}$ , $\Re$ $\partial_z$($\phi$) ($mT$)')
plt.legend(prop={'size': 8})
plt.grid()
#plt.savefig('030ReBZDiff.png', dpi=300, bbox_inches = "tight")
plt.show()

plt.clf()
plt.title(r'Sextupole $\mathrm{\mathbb{Re}}$ Bz and $\partial_z$ ($\mathrm{\mathbb{Re}}$ $\phi$) scalar potential component as a function of z')
plt.plot(zRangeDiff, (np.real(bzDataDiffScale)-np.real(scalarDiff)), label = r'$\mathrm{\mathbb{Re}}$ $\phi$')
plt.xlabel(r'$z$ ($mm$)')
plt.ylabel(r'$\Re$ $B_{z}$ , $\Re$ $\partial_z$($\phi$) ($mT$)')
plt.grid()
#plt.savefig('031ReBZDiffResid.png', dpi=300, bbox_inches = "tight")
plt.show()

plt.clf()
plt.title(r'Sextupole $\mathrm{\mathbb{Re}}$ Bz and $\partial_z$ ($\mathrm{\mathbb{Re}}$ $\phi$) scalar potential component as a function of z')
plt.plot(zRange, np.real(bzData), label = r'$\mathrm{\mathbb{Re}}$ Bz')
plt.plot(zRange, np.real(scalarGrad), label = r'$\partial_z$ ($\mathrm{\mathbb{Re}}$ $\phi$)')
plt.xlabel(r'$z$ ($mm$)')
plt.ylabel(r'$\Re$ $B_{z}$ , $\Re$ $\partial_z$($\phi$) ($mT$)')
plt.legend(prop={'size': 8})
plt.grid()
#plt.savefig('030ReBZDiff.png', dpi=300, bbox_inches = "tight")
plt.show()

plt.clf()
plt.title(r'Sextupole $\mathrm{\mathbb{Re}}$ Bz and $\partial_z$ ($\mathrm{\mathbb{Re}}$ $\phi$) scalar potential component as a function of z')
plt.plot(zRange, (np.real(bzData)-np.real(scalarGrad)), label = r'$\mathrm{\mathbb{Re}}$ $\phi$')
plt.xlabel(r'$z$ ($mm$)')
plt.ylabel(r'$\Re$ $B_{z}$ , $\Re$ $\partial_z$($\phi$) ($mT$)')
plt.grid()
#plt.savefig('031ReBZDiffResid.png', dpi=300, bbox_inches = "tight")
plt.show()
