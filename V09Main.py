###############################################################################
# Version 9.01 Main, Last Modified, 09/12/2021
###############################################################################
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from V09Functions import G, bx, by, bzeta, IntegrateTrapz, NumericalFourierFit
from V09Functions import FFT, FFTFit, ScalarPotential
from V09Params import zLim, zLimPlot
###############################################################################

do_Field_Plots = True

do_Fourier_Analysis = False

do_Fast_Fourier_Transform = False

do_Fourier_Comparisons = False

do_FFT_Scalar_Potential = True

do_2D_Scalar_Potential_Multipole_Plots = True

do_3D_Scalar_Potential_Multipole_Plots = False

###############################################################################

if do_Field_Plots:
    
    figRes = 201
    
    fig1Data = np.zeros([5,201])
    fig1Range = np.linspace(-50,50,201)
    
    for i in range(5):
        for k in range (201):
            fig1Data[i][k] = G((k-100), i)
    
    plt.title('Enge Roll-Off For N-Order Multipoles')
    plt.plot(fig1Range, fig1Data[0,:], label = 'n=0, Dipole')
    plt.plot(fig1Range, fig1Data[1,:], label = 'n=1, Quadrupole')
    plt.plot(fig1Range, fig1Data[2,:], label = 'n=2, Sextupole')
    plt.plot(fig1Range, fig1Data[3,:], label = 'n=3, Octupole')
    plt.plot(fig1Range, fig1Data[4,:], label = 'n=4, Decapole')
    plt.xlabel(r'$z$ ($mm$)')
    plt.ylabel(r'Gradient Function ($mm^{-1}$)')
    plt.grid()
    plt.legend(prop={'size': 8})
    #plt.savefig('001Enge-Roll-Off.png', dpi=300, bbox_inches = "tight")
    plt.show()
    
    Zvalues = np.linspace(-3, 3, 6)
    fig2Range = fig3Range = np.linspace(-0.1, 0.1, figRes)
    fig2Data = np.zeros([6,3,figRes]) 
    fig3Data = np.zeros([6,3,figRes]) 
    
    titlesX = [r'Quadrupole $B_{x}$ Plot',
              r'Sextupole $B_{x}$ Plot',
              r'Octupole $B_{x}$ Plot']
    
    titlesY = [r'Quadrupole $B_{y}$ Plot',
              r'Sextupole $B_{y}$ Plot',
              r'Octupole $B_{y}$ Plot']
    
    savetagX = ['002QuadBXPlot.png',
              '004SextupBXPlot.png',
              '006OctupBXPlot.png']
    
    savetagY = ['003QuadBYPlot.png',
              '005SextupBYPlot.png',
              '007OctupBYPlot.png']
   
    for i in range (3):
        n = i+1
        for j in range(6):
            for k in range (figRes):
                
                fig2Data[j,0,:] = Zvalues[j]
                fig2Data[j,1,k] = fig2Range[k]
                mpc = bx(0.1, fig2Data[j,1,k], fig2Data[j,0,k], n)
                mpc = mpc.real
                fig2Data[j,2,k] = mpc
                
                fig3Data[j,0,:] = Zvalues[j]
                fig3Data[j,1,k] = fig3Range[k]
                mpc = by(fig3Data[j,1,k], 0, fig3Data[j,0,k], n)
                mpc = mpc.real
                fig3Data[j,2,k] = mpc
        
        plt.clf()
        plt.title(titlesX[i])
        plt.plot(fig2Data[0,1,:], fig2Data[0,2,:], label = r'$z$ = -3')
        plt.plot(fig2Data[1,1,:], fig2Data[1,2,:], label = r'$z$ = -1.8')
        plt.plot(fig2Data[2,1,:], fig2Data[2,2,:], label = r'$z$ = -0.6')
        plt.plot(fig2Data[3,1,:], fig2Data[3,2,:], label = r'$z$ = 0.6')
        plt.plot(fig2Data[4,1,:], fig2Data[4,2,:], label = r'$z$ = 1.8')
        plt.plot(fig2Data[5,1,:], fig2Data[5,2,:], label = r'$z$ = 3')
        plt.xlabel(r'$y$ ($mm$)')
        plt.ylabel(r'$B_{x}$ ($mT$)')
        plt.grid()
        plt.legend(prop={'size': 8})
        #plt.savefig(savetagX[i], dpi=300, bbox_inches = "tight")
        plt.show()
        
        plt.clf()
        plt.title(titlesY[i])
        plt.plot(fig3Data[0,1,:], fig3Data[0,2,:], label = r'$z$ = -3')
        plt.plot(fig3Data[1,1,:], fig3Data[1,2,:], label = r'$z$ = -1.8')
        plt.plot(fig3Data[2,1,:], fig3Data[2,2,:], label = r'$z$ = -0.6')
        plt.plot(fig3Data[3,1,:], fig3Data[3,2,:], label = r'$z$ = 0.6')
        plt.plot(fig3Data[4,1,:], fig3Data[4,2,:], label = r'$z$ = 1.8')
        plt.plot(fig3Data[5,1,:], fig3Data[5,2,:], label = r'$z$ = 3')
        plt.xlabel(r'$x$ ($mm$)')
        plt.ylabel(r'$B_{y}$ ($mT$)')
        plt.grid()
        plt.legend(prop={'size': 8})
        #plt.savefig(savetagY[i], dpi=300, bbox_inches = "tight")
        plt.show()

###############################################################################

if do_Fourier_Analysis:
    
    fourierRes = 201
    nFourierCoeffs = 21
    imagpt = 2.3j
    
    hRange, fig3Data = IntegrateTrapz(nFourierCoeffs, fourierRes, 2)
    fig3Data = np.log(np.absolute(fig3Data))
    
    plt.clf()
    plt.title('Fourier Coefficient Values From Numerical Integration')
    plt.plot(hRange, fig3Data, 'ok')
    plt.xlabel(r'$N^{th}$ Harmonic Number')
    plt.ylabel(r'$Ln[Abs(Fourier$ $Coefficient)]$')
    plt.grid()
    #plt.savefig('008NumericalFourierCoeffs.png', dpi=300, bbox_inches = "tight")
    plt.show()
    
    dataGRe = np.zeros([fourierRes])
    dataGIm = np.zeros([fourierRes])
    dataTrapRe = np.zeros([fourierRes])
    dataTrapIm = np.zeros([fourierRes])
    vals = np.linspace(-zLim, zLim, fourierRes)
    
    for i in range(fourierRes):
        
        dataG = G((vals[i]+imagpt), 2)
        dataTrap, zRange = NumericalFourierFit((vals[i]+imagpt), nFourierCoeffs, fourierRes, 2)
        dataGRe[i] = np.real(dataG)
        dataGIm[i] = np.imag(dataG)
        dataTrapRe[i] = np.real(dataTrap)
        dataTrapIm[i] = np.imag(dataTrap) 
    dataDiffRe = dataTrapRe - dataGRe
    dataDiffIm = dataTrapIm - dataGIm
    
    plt.clf()
    plt.title('Numerical Real Fourier Fit Comparison')
    plt.plot(zRange, dataGRe, label = 'G, n = 2')
    plt.plot(zRange, dataTrapRe, label = 'Fourier Fit, n = 2')
    plt.xlabel(r'$z$ ($mm$)')
    plt.ylabel(r'$\Re$ (G) ($mm^{-1}$)')
    plt.legend(prop={'size': 8})
    plt.grid()
    #plt.savefig('009NumericalFourierRe.png', dpi=300, bbox_inches = "tight")
    plt.show()
    
    plt.clf()
    plt.title('Numerical Real Fourier Fit Difference')
    plt.plot(zRange, dataDiffRe)
    plt.xlabel(r'$z$ ($mm$)')
    plt.ylabel(r'Residual of $\Re$ (G) Fit ($mm^{-1}$)')
    plt.grid()
    #plt.savefig('010NumericalFourierReDiff.png', dpi=300, bbox_inches = "tight")
    plt.show()
    
    plt.clf()
    plt.title('Numerical Imaginary Fourier Fit Comparison')
    plt.plot(zRange, dataGIm, label = 'G, n = 2')
    plt.plot(zRange, dataTrapIm, label = 'Fourier Fit, n = 2')
    plt.xlabel(r'$z$ ($mm$)')
    plt.ylabel(r'$\Im$ (G) ($mm^{-1}$)')
    plt.legend(prop={'size': 8})
    plt.grid()
    #plt.savefig('011NumericalFourierIm.png', dpi=300, bbox_inches = "tight")
    plt.show()
    
    plt.clf()
    plt.title('Numerical Imaginary Fourier Fit Difference')
    plt.plot(zRange, dataDiffIm)
    plt.xlabel(r'$z$ ($mm$)')
    plt.ylabel(r'Residual of $\Im$ (G) Fit ($mm^{-1}$)')
    plt.grid()
    #plt.savefig('012NumericalFourierImDiff.png', dpi=300, bbox_inches = "tight")
    plt.show()

###############################################################################

if do_Fast_Fourier_Transform:
    
    fourierRes = 460
    nFourierCoeffs = 45
    imagpt = 2.3j
    
    fig4Data = FFT(nFourierCoeffs, 2)
    fig4Data = np.log(np.absolute(fig4Data))
    fig4Data = fig4Data[1::2]
    fig4Data = fig4Data[0:int((nFourierCoeffs-1)/2)]
    hRange1 = np.linspace(1, int((nFourierCoeffs-2)), int((nFourierCoeffs-1)/2))
    
    plt.clf()
    plt.title('Fourier Coefficient Values From The Fast Fourier Transform')
    plt.plot(hRange1, fig4Data, 'ok')
    plt.xlabel(r'$N^{th}$ Harmonic Number')
    plt.ylabel(r'$Ln[Abs(Fourier$ $Coefficient)]$')
    plt.grid()
    #plt.savefig('013FastFourierCoeffs.png', dpi=300, bbox_inches = "tight")
    plt.show()

    dataGRe2 = np.zeros([fourierRes])
    dataGIm2 = np.zeros([fourierRes])
    dataFFTRe = np.zeros([fourierRes])
    dataFFTIm = np.zeros([fourierRes])
    vals2 = np.linspace(-zLim, zLim, fourierRes)
    
    for i in range(fourierRes):
        
        dataG2 = G((vals2[i]+imagpt), 2)
        dataFFT = FFTFit((vals2[i]+imagpt), nFourierCoeffs, 2)
        dataGRe2[i] = np.real(dataG2)
        dataGIm2[i] = np.imag(dataG2)
        dataFFTRe[i] = np.real(dataFFT)
        dataFFTIm[i] = np.imag(dataFFT) 
    dataDiffRe2 = dataFFTRe - dataGRe2
    dataDiffIm2 = dataFFTIm - dataGIm2
    
    plt.clf()
    plt.title('FFT Real Fourier Fit Comparison')
    plt.plot(vals2, dataGRe2, label = 'G, n = 2')
    plt.plot(vals2, dataFFTRe, label = 'Fast Fourier, n = 2')
    plt.xlabel(r'$z$ ($mm$)')
    plt.ylabel(r'$\Re$ (G) ($mm^{-1}$)')
    plt.legend(prop={'size': 8})
    plt.grid()
    #plt.savefig('014FFTFourierRe.png', dpi=300, bbox_inches = "tight")
    plt.show()
    
    plt.clf()
    plt.title('FFT Real Fourier Fit Difference')
    plt.plot(vals2, dataDiffRe2)
    plt.xlabel(r'$z$ ($mm$)')
    plt.ylabel(r'Residual of $\Re$ (G) Fit ($mm^{-1}$)')
    plt.grid()
    #plt.savefig('015FFTFourierReDiff.png', dpi=300, bbox_inches = "tight")
    plt.show()
    
    plt.clf()
    plt.title('FFT Imaginary Fourier Fit Comparison')
    plt.plot(vals2, dataGIm2, label = 'G, n = 2')
    plt.plot(vals2, dataFFTIm, label = 'Fast Fourier, n = 2')
    plt.xlabel(r'$z$ ($mm$)')
    plt.ylabel(r'$\Im$ (G) ($mm^{-1}$)')
    plt.legend(prop={'size': 8})
    plt.grid()
    #plt.savefig('016FFTFourierIm.png', dpi=300, bbox_inches = "tight")
    plt.show()
    
    plt.clf()
    plt.title('FFT Imaginary Fourier Fit Difference')
    plt.plot(vals2, dataDiffIm2)
    plt.xlabel(r'$z$ ($mm$)')
    plt.ylabel(r'Residual of $\Im$ (G) Fit ($mm^{-1}$)')
    plt.grid()
    #plt.savefig('017FFTFourierImDiff.png', dpi=300, bbox_inches = "tight")
    plt.show()

###############################################################################

if do_Fourier_Comparisons:
    
    plt.clf()
    plt.title('Fourier Coefficients Comparison')
    plt.plot(hRange, fig3Data, 'ok', label = 'Numerical Fourier')
    plt.plot(hRange1, fig4Data, 'or', label = 'Fast Fourier')
    plt.xlabel(r'$N^{th}$ Harmonic Number')
    plt.ylabel(r'$Ln[Abs(Fourier$ $Coefficient)]$')
    plt.legend(prop={'size': 8})
    plt.grid()
    #plt.savefig('018TripFourierCoeffsComparison.png', dpi=300, bbox_inches = "tight")
    plt.show() 
    
    plt.clf()
    plt.title('Numerical Real Fourier Fit Comparison')
    plt.plot(zRange, dataGRe, label = 'G, n = 2')
    plt.plot(zRange, dataTrapRe, label = 'Numerical Fourier, n = 2')
    plt.plot(vals2, dataFFTRe, label = 'Fast Fourier, n = 2')
    plt.xlabel(r'$z$ ($mm$)')
    plt.ylabel(r'$\Re$ (G) ($mm^{-1}$)')
    plt.legend(prop={'size': 8})
    plt.grid()
    #plt.savefig('019TripNumericalFourierRe.png', dpi=300, bbox_inches = "tight")
    plt.show()
    
    plt.clf()
    plt.title('Numerical Real Fourier Fit Difference')
    plt.plot(zRange, dataDiffRe, label = 'Numerical Fourier')
    plt.plot(vals2, dataDiffRe2, label = 'Fast Fourier')
    plt.xlabel(r'$z$ ($mm$)')
    plt.ylabel(r'Residual of $\Re$ (G) Fit ($mm^{-1}$)')
    plt.legend(prop={'size': 8})
    plt.grid()
    #plt.savefig('020TripNumericalFourierReDiff.png', dpi=300, bbox_inches = "tight")
    plt.show()
    
    plt.clf()
    plt.title('Numerical Imaginary Fourier Fit Comparison')
    plt.plot(zRange, dataGIm, label = 'G, n = 2')
    plt.plot(zRange, dataTrapIm, label = 'Numerical Fourier, n = 2')
    plt.plot(vals2, dataFFTIm, label = 'Fast Fourier, n = 2')
    plt.xlabel(r'$z$ ($mm$)')
    plt.ylabel(r'$\Im$ (G) ($mm^{-1}$)')
    plt.legend(prop={'size': 8})
    plt.grid()
    #plt.savefig('021TripNumericalFourierIm.png', dpi=300, bbox_inches = "tight")
    plt.show()
    
    plt.clf()
    plt.title('Numerical Imaginary Fourier Fit Difference')
    plt.plot(zRange, dataDiffIm, label = 'Numerical Fourier')
    plt.plot(vals2, dataDiffIm2, label = 'Fast Fourier')
    plt.xlabel(r'$z$ ($mm$)')
    plt.ylabel(r'Residual of $\Im$ (G) Fit ($mm^{-1}$)')
    plt.legend(prop={'size': 8})
    plt.grid()
    #plt.savefig('022TripNumericalFourierImDiff.png', dpi=300, bbox_inches = "tight")
    plt.show()

###############################################################################

if do_FFT_Scalar_Potential:
    
    fourierRes = 1001
    nFourierCoeffs = 45
    
    BXfactor = 4.05
    BYfactor = 9.7
    diffFactor = 177

    x = 0.5
    y = 0.1
    zRange = np.linspace(-zLimPlot, zLimPlot, fourierRes)
    zRangeDiff = np.linspace((-zLimPlot+(1/(2*fourierRes))),(zLimPlot-(1/(2*fourierRes))),(fourierRes-1)) # generate the new z scale.
    
    n = 2 # Multipole Number - Sextupole
    
    bxData = np.zeros([fourierRes], dtype = np.complex_)
    byData = np.zeros([fourierRes], dtype = np.complex_)
    bzData = np.zeros([fourierRes], dtype = np.complex_)
    bzDataDiffScale = np.zeros([fourierRes-1], dtype = np.complex_)
    scalarPotData = np.zeros([fourierRes], dtype = np.complex_)
    
    for i in range(fourierRes):
        
        bxDataPoint = bx(x, y, zRange[i], n)
        byDataPoint = by(x, y, zRange[i], n)
        bzDataPoint = bzeta(x, y, zRange[i], n)
        
        bxData[i] = bxDataPoint
        byData[i] = byDataPoint
        bzData[i] = bzDataPoint
       
    for i in range(fourierRes-1):
        
        bzDataPointDiff = bzeta(x, y, zRangeDiff[i], n)
        bzDataDiffScale[i] = bzDataPointDiff
    
    scalarPotData = ScalarPotential(x, y, zRange, nFourierCoeffs, n) # Compute first for FourierRes
    scalarDiff = np.diff(scalarPotData, n=1)
    scalarDiff = scalarDiff*diffFactor # Differentiate numerically  which alters the z scale.
    
    plt.clf()
    plt.title(r'Sextupole $\mathrm{\mathbb{Re}}$ Bz Component as a Function of z')
    plt.plot(zRange, np.real(bzData), label = r'$\mathrm{\mathbb{Re}}$ Bz')
    plt.xlabel(r'$z$ ($mm$)')
    plt.ylabel(r'$\Re$ $B_{z}$ ($mT$)')
    plt.grid()
    #plt.savefig('023Bz.png', dpi=300, bbox_inches = "tight")
    plt.show()
    
    plt.clf()
    plt.title(r'Sextupole $\mathrm{\mathbb{Re}}$ Scalar Potential as a Function of z')
    plt.plot(zRange, np.real(scalarPotData), label = r'$\mathrm{\mathbb{Re}}$ $\phi$')
    plt.xlabel(r'$z$ ($mm$)')
    plt.ylabel(r'$\Re$ $\phi$ ($mT$ $mm$)')
    plt.grid()
    #plt.savefig('024ReScalarPotential.png', dpi=300, bbox_inches = "tight")
    plt.show()
    
    plt.clf()
    plt.title(r'Sextupole $\mathrm{\mathbb{Re}}$ $\partial_z$ scalar potential component as a function of z')
    plt.plot(zRangeDiff, np.real(scalarDiff), label = r'$\partial_z$ ($\mathrm{\mathbb{Re}}$ $\phi$)')
    plt.xlabel(r'$z$ ($mm$)')
    plt.ylabel(r'$\Re$ $\partial_z$($\phi$) ($mT$)')
    plt.grid()
    #plt.savefig('025ScalarPotDiff.png', dpi=300, bbox_inches = "tight")
    plt.show()
    
    plt.clf()
    plt.title(r'Sextupole $\mathrm{\mathbb{Re}}$ Bx and $\partial_x$ ($\mathrm{\mathbb{Re}}$ $\phi$) component as a function of z')
    plt.plot(zRange, np.real(bxData), label = r'$\mathrm{\mathbb{Re}}$ Bx')
    plt.plot(zRange, (BXfactor*np.real(scalarPotData)), label = r'$\mathrm{\mathbb{Re}}$ $\phi$')
    plt.xlabel(r'$z$ ($mm$)')
    plt.ylabel(r'$\Re$ $B_{x}$ , $\Re$ $\partial_x$($\phi$) ($mT$)')
    plt.legend(prop={'size': 8})
    plt.grid()
    #plt.savefig('026ReBXDiff.png', dpi=300, bbox_inches = "tight")
    plt.show()
    
    plt.clf()
    plt.title(r'Sextupole $\mathrm{\mathbb{Re}}$ Bx and $\partial_x$ ($\mathrm{\mathbb{Re}}$ $\phi$) component as a function of z')
    plt.plot(zRange, (np.real(bxData)-(BXfactor*np.real(scalarPotData))), label = r'$\mathrm{\mathbb{Re}}$ Bx')
    plt.xlabel(r'$z$ ($mm$)')
    plt.ylabel(r'$\Re$ $B_{x}$ , $\Re$ $\partial_x$($\phi$) ($mT$)')
    plt.grid()
    #plt.savefig('027ReBXDiffResid.png', dpi=300, bbox_inches = "tight")
    plt.show()
    
    plt.clf()
    plt.title(r'Sextupole $\mathrm{\mathbb{Re}}$ By and $\partial_y$ ($\mathrm{\mathbb{Re}}$ $\phi$) component as a function of z')
    plt.plot(zRange, np.real(byData), label = r'$\mathrm{\mathbb{Re}}$ By')
    plt.plot(zRange, (BYfactor*np.real(scalarPotData)), label = r'$\mathrm{\mathbb{Re}}$ $\phi$')
    plt.xlabel(r'$z$ ($mm$)')
    plt.ylabel(r'$\Re$ $B_{y}$ , $\Re$ $\partial_y$($\phi$) ($mT$)')
    plt.legend(prop={'size': 8})
    plt.grid()
    #plt.savefig('028ReBYDiff.png', dpi=300, bbox_inches = "tight")
    plt.show()
    
    plt.clf()
    plt.title(r'Sextupole $\mathrm{\mathbb{Re}}$ By and $\partial_y$ ($\mathrm{\mathbb{Re}}$ $\phi$) component as a function of z')
    plt.plot(zRange, (np.real(byData)-(BYfactor*np.real(scalarPotData))), label = r'$\mathrm{\mathbb{Re}}$ By')
    plt.xlabel(r'$z$ ($mm$)')
    plt.ylabel(r'$\Re$ $B_{y}$ , $\Re$ $\partial_y$($\phi$) ($mT$)')
    plt.grid()
    #plt.savefig('029ReBYDiffResid.png', dpi=300, bbox_inches = "tight")
    plt.show()
    
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

###############################################################################

if do_2D_Scalar_Potential_Multipole_Plots:
    
    fourierRes = 51 #201
    nFourierCoeffs = 35 #45

    
    
    xVals = np.linspace(-1.5, 1.5, fourierRes)
    yVals = np.linspace(-1.5, 1.5, fourierRes)
    X , Y = np.meshgrid(xVals,yVals)
    datasetDipole = np.zeros([fourierRes,fourierRes])
    datasetQuadrupole = np.zeros([fourierRes,fourierRes])
    datasetSextupole = np.zeros([fourierRes,fourierRes])
    datasetOctupole = np.zeros([fourierRes,fourierRes])
    
    for k in range(fourierRes):
        print('Run ', k ,' out of ', fourierRes)
        for i in range(fourierRes):
            
            x = X[0,k]
            y = Y[i,0]
            
            dataPoint0 = ScalarPotential(x, y, 0, nFourierCoeffs, 0)
            dataPoint1 = ScalarPotential(x, y, 0, nFourierCoeffs, 1)
            dataPoint2 = ScalarPotential(x, y, 0, nFourierCoeffs, 2)
            dataPoint3 = ScalarPotential(x, y, 0, nFourierCoeffs, 3)
            
            datasetDipole[k,i] = dataPoint0
            datasetQuadrupole[k,i] = dataPoint1
            datasetSextupole[k,i] = dataPoint2
            datasetOctupole[k,i] = dataPoint3
    
    plt.clf()
    plt.pcolormesh(X, Y, datasetDipole, cmap = cm.seismic) 
    plt.title('Dipole x,y Scalar Potential Heat Map, z = 0')
    plt.xlabel('x (mm)')
    plt.ylabel('y (mm)')
    #plt.savefig('032dipole.png', dpi=300, bbox_inches = "tight")
    plt.show()  
    
    plt.clf()
    plt.pcolormesh(X, Y, datasetQuadrupole, cmap = cm.seismic) 
    plt.title('Quadrupole x,y Scalar Potential Heat Map, z = 0')
    plt.xlabel('x (mm)')
    plt.ylabel('y (mm)')
    #plt.savefig('033quadruople.png', dpi=300, bbox_inches = "tight")
    plt.show()
    
    plt.clf()
    plt.pcolormesh(X, Y, datasetSextupole, cmap = cm.seismic) 
    plt.title('Sextupole x,y Scalar Potential Heat Map, z = 0')
    plt.xlabel('x (mm)')
    plt.ylabel('y (mm)')
    #plt.savefig('034sextupole.png', dpi=300, bbox_inches = "tight")
    plt.show()
    
    plt.clf()
    plt.pcolormesh(X, Y, datasetOctupole, cmap = cm.seismic) 
    plt.title('Octupole x,y Scalar Potential Heat Map, z = 0')
    plt.xlabel('x (mm)')
    plt.ylabel('y (mm)')
    #plt.savefig('035octupole.png', dpi=300, bbox_inches = "tight")
    plt.show()
    
    plt.clf()
    plt.contour(X, Y, datasetQuadrupole, levels = [-0.5, -0.4, -0.3, -0.25, -0.2, 0.2, 0.25, 0.3, 0.4, 0.5])#, cmap = cm.seismic) 
    plt.title('Quadrupole x,y Scalar Potential Heat Map, z = 0')
    plt.xlabel('x (mm)')
    plt.ylabel('y (mm)')
    #plt.savefig('036ContourQuadruople.png', dpi=300, bbox_inches = "tight")
    plt.show()
    
    plt.clf()
    plt.contour(X, Y, datasetSextupole, levels = [-0.5, -0.4, -0.3, -0.25, -0.2, 0.2, 0.25, 0.3, 0.4, 0.5])#, cmap = cm.seismic) 
    plt.title('Sextupole x,y Scalar Potential Heat Map, z = 0')
    plt.xlabel('x (mm)')
    plt.ylabel('y (mm)')
    #plt.savefig('037ContourSextupole.png', dpi=300, bbox_inches = "tight")
    plt.show()
    
    plt.clf()
    plt.contour(X, Y, datasetOctupole, levels = [-0.4, -0.3, -0.25, -0.2, -0.15, 0.15, 0.2, 0.25, 0.3, 0.4])#, cmap = cm.seismic) 
    plt.title('Octupole x,y Scalar Potential Heat Map, z = 0')
    plt.xlabel('x (mm)')
    plt.ylabel('y (mm)')
    #plt.savefig('038ContourOctupole.png', dpi=300, bbox_inches = "tight")
    plt.show()
    
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    surf = ax.plot_surface(X, Y, np.absolute(datasetQuadrupole), cmap=cm.coolwarm,
                           linewidth=0, antialiased=False)
    plt.savefig('039_3DContourOctupole.png', dpi=300, bbox_inches = "tight")
    plt.show()

    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    surf = ax.plot_surface(X, Y, np.absolute(datasetSextupole), cmap=cm.coolwarm,
                           linewidth=0, antialiased=False)
    plt.savefig('040_3DContourOctupole.png', dpi=300, bbox_inches = "tight")
    plt.show()

    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    surf = ax.plot_surface(X, Y, np.absolute(datasetOctupole), cmap=cm.coolwarm,
                           linewidth=0, antialiased=False)
    plt.savefig('041_3DContourOctupole.png', dpi=300, bbox_inches = "tight")
    plt.show()

###############################################################################

if do_3D_Scalar_Potential_Multipole_Plots:
    
    fourierRes = 25#51 #201
    nFourierCoeffs = 19#35 #45

    
    
    xVals = np.linspace(-1.5, 1.5, fourierRes)
    yVals = np.linspace(-1.5, 1.5, fourierRes)
    zVals = np.linspace(-0.5, 0.5, 3)
    X , Y = np.meshgrid(xVals,yVals)
    datasetDipole = np.zeros([3,fourierRes,fourierRes])
    datasetQuadrupole = np.zeros([3,fourierRes,fourierRes])
    datasetSextupole = np.zeros([3,fourierRes,fourierRes])
    datasetOctupole = np.zeros([3,fourierRes,fourierRes])
    
    for m in range (3):
        for k in range(fourierRes):
            print('Run ', k ,' out of ', fourierRes)
            for i in range(fourierRes):
                
                x = X[0,k]
                y = Y[i,0]
                
                #dataPoint0 = ScalarPotential(x, y, zVals[m], nFourierCoeffs, 0)
                #dataPoint1 = ScalarPotential(x, y, zVals[m], nFourierCoeffs, 1)
                dataPoint2 = ScalarPotential(x, y, zVals[m], nFourierCoeffs, 2)
                #dataPoint3 = ScalarPotential(x, y, zVals[m], nFourierCoeffs, 3)
                
                #datasetDipole[m,k,i] = dataPoint0
                #datasetQuadrupole[m,k,i] = dataPoint1
                datasetSextupole[m,k,i] = dataPoint2
                #datasetOctupole[m,k,i] = dataPoint3


    
    zz = np.linspace(-1, 1, 3)
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    surf = ax.plot_surface(X, np.absolute(datasetSextupole[1,:,:]), zz, cmap=cm.coolwarm,
                           linewidth=0, antialiased=False)
    plt.show()
    plt.clf()
    '''
    plt.plot_surface(X, Y, np.absolute(datasetSextupole[1,:,:]), cmap=cm.coolwarm,
                    linewidth=0, antialiased=False)
    plt.show()
    '''