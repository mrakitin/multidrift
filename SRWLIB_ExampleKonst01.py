#############################################################################
# SRWLIB Example for Kon's R&D #1: Simulating propagation of a Gaussian X-ray beam through a Beamline containing a ZP and Other Stuff
# v 0.01
#############################################################################

from __future__ import print_function #Python 2.7 compatibility
from srwlib import *
from uti_plot import * #required for plotting
#import os
#import sys
import time

print('SRWLIB Example for Kon\'s R&D #1')
print('Simulating propagation of a Gaussian X-ray beam through a Beamline containing a ZP and Other Stuff')

#**********************Input Parameters and Structures
#***********Folder and Data File Names
strDataFolderName = 'data_konst' #data sub-folder name
strIntInitOutFileName01 = 'res_int_in.dat' #initial wavefront intensity distribution output file name
strIntPropOutFileName01 = 'res_int_prop.dat' #propagated wavefront intensity distribution output file name

#***********Gaussian Beam Source
GsnBm = SRWLGsnBm() #Gaussian Beam structure (just parameters)
GsnBm.x = 0 #Transverse Positions of Gaussian Beam Center at Waist [m]
GsnBm.y = 0
GsnBm.z = 0 #Longitudinal Position of Waist [m]
GsnBm.xp = 0 #Average Angles of Gaussian Beam at Waist [rad]
GsnBm.yp = 0
GsnBm.avgPhotEn = 930 #Photon Energy [eV]
GsnBm.pulseEn = 0.001 #Energy per Pulse [J] - to be corrected
GsnBm.repRate = 1 #Rep. Rate [Hz] - to be corrected
GsnBm.polar = 1 #1- linear hoirizontal
GsnBm.sigX = 2.e-09/2.35 #23e-06/2.35 #Horiz. RMS size at Waist [m]
GsnBm.sigY = GsnBm.sigX #Vert. RMS size at Waist [m]

constConvRad = 1.23984186e-06/(4*3.1415926536)
rmsAngDiv = constConvRad/(GsnBm.avgPhotEn*GsnBm.sigX) #RMS angular divergence [rad]
print('RMS Source Size:', round(GsnBm.sigX*1.e+06, 3), 'microns; RMS Divergence:', round(rmsAngDiv*1.e+03, 3), 'mrad')

GsnBm.sigT = 10e-15 #Pulse duration [fs] (not used?)
GsnBm.mx = 0 #Transverse Gauss-Hermite Mode Orders
GsnBm.my = 0

#***********Initial Wavefront
wfr = SRWLWfr() #Initial Electric Field Wavefront
wfr.allocate(1, 100, 100) #Numbers of points vs Photon Energy (1), Horizontal and Vertical Positions (dummy)
wfr.mesh.zStart = 0 #Longitudinal Position [m] at which initial Electric Field has to be calculated, i.e. the position of the first optical element
#The above parameter will be re-defined below.

wfr.mesh.eStart = GsnBm.avgPhotEn #Initial Photon Energy [eV]
wfr.mesh.eFin = GsnBm.avgPhotEn #Final Photon Energy [eV]

wfr.unitElFld = 2 #Electric field units: 0- arbitrary, 1- sqrt(Phot/s/0.1%bw/mm^2), 2- sqrt(J/eV/mm^2) or sqrt(W/mm^2), depending on representation (freq. or time)

distSrc_M1 = wfr.mesh.zStart - GsnBm.z
#Horizontal and Vertical Position Range for the Initial Wavefront calculation
#can be used to simulate the First Aperture (of M1)
firstHorAp = 0.5e-03 #2.e-03 #8.*rmsAngDiv*distSrc_M1 #[m]
firstVertAp = firstHorAp #[m]

wfr.mesh.xStart = -0.5*firstHorAp #Initial Horizontal Position [m]
wfr.mesh.xFin = 0.5*firstHorAp #Final Horizontal Position [m]
wfr.mesh.yStart = -0.5*firstVertAp #Initial Vertical Position [m]
wfr.mesh.yFin = 0.5*firstVertAp #Final Vertical Position [m]

sampFactNxNyForProp = 1.5 #2.5 #1.25 #2 #sampling factor for adjusting nx, ny (effective if > 0)
arPrecPar = [sampFactNxNyForProp]

wfr.partBeam.partStatMom1.x = GsnBm.x #Some information about the source in the Wavefront structure
wfr.partBeam.partStatMom1.y = GsnBm.y
wfr.partBeam.partStatMom1.z = GsnBm.z
wfr.partBeam.partStatMom1.xp = GsnBm.xp
wfr.partBeam.partStatMom1.yp = GsnBm.yp

#***********Optical Elements and Propagation Parameters
#Sequence of Optical Elements:
#   <Aperture of ZP>
#   <Zone Plate (ZP)>
#   <Drift from ZP to Slit>
#   <Slit>
#   <Drift from Slit to ...>

#ZP Parameters
cenPhotEnZP = GsnBm.avgPhotEn #ZP central Photon Energy
diamZP = firstHorAp #ZP Diameter
numZonesZP = 1000 #ZP Number of Zones
focLenZP = 2.0164e+5*diamZP*diamZP*cenPhotEnZP/numZonesZP #Estimated Focal Length
print('   ZP Focal Length:', focLenZP, 'm')

distZP_S = 2.5 #Distance from ZP to Slit

distSrc_ZP = distZP_S*focLenZP/(distZP_S - focLenZP) #Distance from ZP to Slit
print('   ZP position resp. to \'Source\':', distSrc_ZP, 'm')

wfr.mesh.zStart = distSrc_ZP

sizeSx = 1.e-06 #Slit sizes
sizeSy = sizeSx
distS_Obs = 0.05 #Distance from Slit to Observation plane

#Aperture of M1
opAZP = SRWLOptA('c', 'a', diamZP, diamZP)

#Zone Plate
opZP = SRWLOptZP(_nZones=1000, _rn=0.5*diamZP, _thick=50e-06, _delta1=1e-06, _atLen1=0.1, _delta2=0, _atLen2=1e-06)

#Drift from Zp to Slit
opZP_S = SRWLOptD(distZP_S) 

#Slit
opS = SRWLOptA('r', 'a', sizeSx, sizeSy)

#Drift after Slit
opS_O = SRWLOptD(distS_Obs)

#Wavefront Propagation Parameters:
#[0]: Auto-Resize (1) or not (0) Before propagation
#[1]: Auto-Resize (1) or not (0) After propagation
#[2]: Relative Precision for propagation with Auto-Resizing (1. is nominal)
#[3]: Allow (1) or not (0) for semi-analytical treatment of the quadratic (leading) phase terms at the propagation
#[4]: Do any Resizing on Fourier side, using FFT, (1) or not (0)
#[5]: Horizontal Range modification factor at Resizing (1. means no modification)
#[6]: Horizontal Resolution modification factor at Resizing
#[7]: Vertical Range modification factor at Resizing
#[8]: Vertical Resolution modification factor at Resizing
#[9]: Type of wavefront Shift before Resizing (not yet implemented)
#[10]: New Horizontal wavefront Center position after Shift (not yet implemented)
#[11]: New Vertical wavefront Center position after Shift (not yet implemented)
#        [0][1][2] [3][4] [5] [6]  [7]  [8] [9][10][11] 
ppAZP =  [0, 0, 1., 0, 0, 1.,  1.,  1.,  1.,  0, 0, 0]
ppZP =   [0, 0, 1., 0, 0, 1.,  1.,  1.,  1.,  0, 0, 0]
ppZP_S = [0, 0, 1., 1, 0, 1.3, 1.,  1.3, 1.,  0, 0, 0]
ppS =    [0, 0, 1., 0, 0, 1.,  1.,  1.,  1.,  0, 0, 0]
ppS_O =  [0, 0, 1., 1, 0, 1.,  1.,  1.,  1.,  0, 0, 0]
ppFin =  [0, 0, 1., 0, 1, 0.05,3.,  0.05,3.,  0, 0, 0] #Post-propagation resizing parameters

#"Beamline" - a sequenced Container of Optical Elements (together with the corresponding wavefront propagation parameters,
#and the "post-propagation" wavefront resizing parameters for better viewing).
opBL = SRWLOptC([opAZP, opZP, opZP_S],# opS, opS_O],
                [ppAZP, ppZP, ppZP_S, ppFin])#, ppS, ppS_O, ppFin])

#**********************Calculation
#***********Initial Wavefront of Gaussian Beam

print('Calculating initial wavefront and eventually extracting its intensity...', end='')
t0 = time.time();
srwl.CalcElecFieldGaussian(wfr, GsnBm, arPrecPar)
print('done in', round(time.time() - t0), 's')

print('   Numbers of points in the initial wavefront: nx=', wfr.mesh.nx, ' ny=', wfr.mesh.ny)

print('Extracting the initial wavefront intensity and eventually saving it to a file...', end='')
t0 = time.time();
mesh0 = deepcopy(wfr.mesh)
arI0 = array('f', [0]*mesh0.nx*mesh0.ny) #"flat" array to take 2D intensity data
srwl.CalcIntFromElecField(arI0, wfr, 6, 0, 3, mesh0.eStart, 0, 0) #extracts intensity
#srwl_uti_save_intens_ascii(arI0, mesh0, os.path.join(os.getcwd(), strDataFolderName, strIntInitOutFileName01), 0,
#                           ['Photon Energy', 'Horizontal Position', 'Vertical Position', 'Spectral Fluence'], _arUnits=['eV', 'm', 'm', 'J/eV/mm^2'])
print('done in', round(time.time() - t0), 's')

#***********Wavefront Propagation
print('Propagating wavefront...', end='')
t0 = time.time();
srwl.PropagElecField(wfr, opBL)
print('done in', round(time.time() - t0), 's')

print('Extracting the propagated wavefront intensity and saving it to a file...', end='')
mesh1 = deepcopy(wfr.mesh)
arI1 = array('f', [0]*mesh1.nx*mesh1.ny) #"flat" array to take 2D intensity data
srwl.CalcIntFromElecField(arI1, wfr, 6, 0, 3, mesh1.eStart, 0, 0) #extracts intensity
srwl_uti_save_intens_ascii(arI1, mesh1, os.path.join(os.getcwd(), strDataFolderName, strIntPropOutFileName01), 0,
                           ['Photon Energy', 'Horizontal Position', 'Vertical Position', 'Spectral Fluence'], _arUnits=['eV', 'm', 'm', 'J/eV/mm^2'])
print('done')

#**********************Plotting Results (requires 3rd party graphics package)
print('   Plotting the results (blocks script execution; close any graph windows to proceed) ... ', end='')
plotMesh0x = [mesh0.xStart, mesh0.xFin, mesh0.nx]
plotMesh0y = [mesh0.yStart, mesh0.yFin, mesh0.ny]
uti_plot2d1d(arI0, plotMesh0x, plotMesh0y, labels=['Horizontal Position', 'Vertical Position', 'Intensity Before Propagation'], units=['m', 'm', 'J/eV/mm^2'])

plotMesh1x = [mesh1.xStart, mesh1.xFin, mesh1.nx]
plotMesh1y = [mesh1.yStart, mesh1.yFin, mesh1.ny]
uti_plot2d1d(arI1, plotMesh1x, plotMesh1y, labels=['Horizontal Position', 'Vertical Position', 'Intensity After Propagation'], units=['m', 'm', 'J/eV/mm^2'])

uti_plot_show() #show all graphs (blocks script execution; close all graph windows to proceed)
print('done')

