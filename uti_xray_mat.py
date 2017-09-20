# -*- coding: utf-8 -*-
"""
Created on Fri Apr 05 16:19:42 2013
package contains standard x-ray related functions from xfuncs package (version 2.0 from 03/27/2015).
changed all defaults to SI units!
Re-packaged for use in SRW 05/29/2015
Needs to be linked to database!!!
@author: Lutz Wiegart, O.Chubar
"""
#****************************************************************************

from __future__ import print_function #Python 2.7 compatibility
#import pylab as pl
import uti_io #OC150715
import uti_math
#import numpy as np
#from os import listdir
#from os.path import isfile, join
import os #OC150715
import re
from array import *

# path to X-ray database files
#datapath='E:/lutz/srw_python_virt_bl_chx/Xray_database_SRW/'
datapath = os.path.join(os.getcwd(), 'data_xray_mat') #OC150715

#****************************************************************************
def get_Lambda(E,u='SI'):
    """
    calculates X-ray wavelength as a function of Energy [eV] in optional units.
    Syntax: getLambda(E,u), 
    where E=X-ray energy; optional: u= 'A','nm','um','cm','mm','m','SI' (='m'), default in the absence of u: 'SI'  
    
    """
    hPlank=6.62606876e-34;
    cvac=2.99792458e8;
    Qelectron=1.602176463e-19;
    scale=1
    l=hPlank*cvac/(E*Qelectron);
    if u is 'A':
        scale=1e10;return l*scale # Angstroem
    elif u is 'nm':
        scale=1e9; return l*scale # nm
    elif u is 'um':
        scale=1e6; return l*scale # um
    elif u is 'mm':
        scale=1e3; return l*scale # mm
    elif u is 'cm':
        scale=1e2; return l*scale # cm
    elif u is 'm' or u is 'SI':
        scale=1; return l*scale
    else:
        #print 'invalid option, type "get_Lambda(\'?\')" for available options and syntax'
        print('invalid option, type "get_Lambda(\'?\')" for available options and syntax')

#****************************************************************************
def get_ac(material,E=8000):
    """
    by LW 10/03/2010
    function calculates the critical angle for total external reflection as a function of
    the material and the X-ray energy according to ac=sqrt(2*delta)
    index of refraction is a .dat file from http://henke.lbl.gov/optical_constants/getdb2.html
    (energy range: 2-30keV,delete the header % lines, name the file n_material.dat) % 
    calling sequence: ac=get_ac(material,E) where ac: critial angle in degrees, E [eV] (default: 8000eV)
    type get_ac(\'materilal?\') to show list of supported materials"
    """
    
    #get list_of supported materials from data file directory:
    #xdatafiles = [ f for f in listdir(datapath) if isfile(join(datapath,f)) ]
    xdatafiles = [ f for f in os.listdir(datapath) if os.path.isfile(os.path.join(datapath,f)) ] #OC150715
    
    name=[]
    #for i in range(0, np.size(xdatafiles)):
    for i in range(0, len(xdatafiles)): #OC150715
        m=re.search('(?<=n_)\w+', xdatafiles[i])
        if m is not None:
            name.append(m.group(0))             

    #E=np.array(E)
    if material in name:
        #loadn=datapath+'n_'+material+'.dat'
        #n=pl.loadtxt(loadn,comments='%')
        loadn = os.path.join(datapath, 'n_'+material+'.dat') #OC150715
        n = uti_io.read_ascii_data_cols(loadn, ' ', _i_col_start=0, _i_col_end=1, _n_line_skip=2) #OC150715

        #if np.min(E)>=np.min(n[:,0]) and np.max(E)<=np.max(n[:,0]):
        #    d=np.interp(E,n[:,0],n[:,1])
        #    return np.degrees(np.sqrt(2*d))
        #else: raise xfuncs_Exception ('error: energy '+"%3.4f" %E +'[eV] out of range ('+"%3.4f" % np.min(n[:,0])+'=<E<='+"%3.4f" % np.max(n[:,0])+'eV)')

        arE = n[0]; arDelta = n[1]
        ####################################################
            

    elif material=='material?':
        #print 'list of supported materials (based on data files in directory '+datapath+':'
        #print name
        print('list of supported materials (based on data files in directory '+datapath+':')
        print(name)
        
    else: raise xfuncs_Exception ('error: non recognized material, please create index of refraction file first. Type "get_ac?" for instructions; type get_ac("material?") for list of supported materials' )

#****************************************************************************      
#def get_n(material, E=8000):
def get_refr(material, E=8000):
    """
    "by LW 07/04/2011 function get the complex index of refraction from stored data file,
    index of refraction is a .dat file from http://henke.lbl.gov/optical_constants/getdb2.html 
    (energy range: 2-30keV,delete the header lines, name the file n_material.dat) 
    calling sequence: n=get_n(material,E) where n is the complex refractive index detlta-i*beta, E: X-ray energy in eV, default: 8000eV"
    """
    #get list_of supported materials from data file directory:
    #xdatafiles = [ f for f in listdir(datapath) if isfile(join(datapath,f)) ]
    xdatafiles = [ f for f in os.listdir(datapath) if os.path.isfile(os.path.join(datapath,f)) ] #OC150715
    
    name=[]
    #for i in range(0, np.size(xdatafiles)):
    for i in range(0, len(xdatafiles)): #OC160715
        m=re.search('(?<=n_)\w+', xdatafiles[i])
        if m is not None:
            name.append(m.group(0))             
    
    #E=np.array(E)
    if material in name:
        #loadn=datapath+'n_'+material+'.dat'
        #n=pl.loadtxt(loadn,comments='%')
        loadn = os.path.join(datapath, 'n_'+material+'.dat') #OC160715
        n = uti_io.read_ascii_data_cols(loadn, ' ', _i_col_start=0, _i_col_end=2, _n_line_skip=2) #OC160715
        
        #if np.min(E)>=np.min(n[:,0]) and np.max(E)<=np.max(n[:,0]):
        #    d=np.interp(E,n[:,0],n[:,1])
        #    b=np.interp(E,n[:,0],n[:,2])
        #    return d-1j*b
        #else: raise xfuncs_Exception ('error: energy '+"%3.4f" %E +'[keV] out of range ('+"%3.4f" % np.min(n[:,0])+'=<E<='+"%3.4f" % np.max(n[:,0])+'eV)')

        arE = n[0]; arDelta = n[1]; arBeta = n[2]
        d = uti_math.interp_1d_var(E, arE, arDelta)
        b = uti_math.interp_1d_var(E, arE, arBeta)
        
        nE_mi_1 = len(arE) - 1
        if((E < arE[0]) or (E > arE[nE_mi_1])):
            print('Warning: energy '+"%3.4f" %E +'[keV] out of range ('+"%3.4f" % arE[0]+'=<E<='+"%3.4f" % arE[nE_mi_1]+'eV)')

        return d - 1j*b #Check the sign of imaginary part!

    elif material=='material?':
        #print 'list of supported materials (based on data files in directory '+datapath+':'
        #print name
        print('list of supported materials (based on data files in directory '+datapath+':')
        print(name)
    #else: raise xfuncs_Exception ( 'error: non recognized material, please create index of refraction file first. Type "get_n?" for instructions; type get_n("material?") for list of supported materials')
    else: raise uti_xray_mat_Exception('error: non recognized material, please create index of refraction file first. Type "get_n?" for instructions; type get_n("material?") for list of supported materials')

#****************************************************************************
#def get_mu(material, E=8000):
def get_att_len(material, E=8000):
    """
    by LW 07/04/2011
    function gets the attenuation length from stored data file, 
    attenuation length is a .dat file from http://henke.lbl.gov/optical_constants/getdb2.html
    (energy range: 2-30keV,delete the header lines or comment with '%', name the file n_material.dat)
    calling sequence: mu=get_mu(material,E) where mu [meter] is the 1/e attenuation length, E: X-ray energy in eV, default: 8000eV'
    """
    #get list_of supported materials from data file directory:
    #xdatafiles = [ f for f in listdir(datapath) if isfile(join(datapath,f)) ]
    xdatafiles = [ f for f in os.listdir(datapath) if os.path.isfile(os.path.join(datapath,f)) ] #OC150715
    
    name=[]
    #for i in range(0, np.size(xdatafiles)):
    for i in range(0, len(xdatafiles)): #OC100915
        mm=re.search('(?<=mu_)\w+', xdatafiles[i])
        if mm is not None:
            name.append(mm.group(0))             
    
    #E=np.array(E)
    if material in name:
        #loadn=datapath+'mu_'+material+'.dat'
        loadm = os.path.join(datapath, 'mu_'+material+'.dat') #OC100915
        #m=pl.loadtxt(loadn,comments='%')
        m = uti_io.read_ascii_data_cols(loadm, ' ', _i_col_start=0, _i_col_end=1, _n_line_skip=2) #OC100915
        
        #if np.min(E)>=np.min(m[:,0]) and np.max(E)<=np.max(m[:,0]):
        #    mu=np.interp(E,m[:,0],m[:,1])
        #    return mu*1e-6
        #else: raise xfuncs_Exception ('error: energy '+"%3.4f" %E +'[keV] out of range ('+"%3.4f" % np.min(m[:,0])+'=<E<='+"%3.4f" % np.max(m[:,0])+'eV)')

        arE = m[0]; arAttLen = m[1]
        attLen = (1e-6)*uti_math.interp_1d_var(E, arE, arAttLen)

        nE_mi_1 = len(arE) - 1
        if((E < arE[0]) or (E > arE[nE_mi_1])):
            print('Warning: energy '+"%3.4f" %E +'[keV] out of range ('+"%3.4f" % arE[0]+'=<E<='+"%3.4f" % arE[nE_mi_1]+'eV)')

        return attLen

    elif material=='material?':
        #print 'list of supported materials (based on data files in directory '+datapath+':'
        #print name
        print('list of supported materials (based on data files in directory '+datapath+':')
        print(name)
    else: raise xfuncs_Exception ('error: non recognized material, please create index of refraction file first. Type get_mu("?") for instructions; type get_n("material?") for list of supported materials')

#****************************************************************************
def get_T(material,E=8000,l=1):
    """
    by LW 10/03/2010, 
    function calculates the transmission as a function of the material and the X-ray energy according to e^(-mul),
    where mu=4pi/lambda*beta 
    index of refraction is a .dat file from http://henke.lbl.gov/optical_constants/getdb2.html 
    (energy range: 2-30keV,delete the header lines, name the file n_material.dat) 
    calling sequence: T=get_T(material,E,l) 
    where T: transmission, material: E: X-ray energy in eV (default: 8000eV), l: thickness of the material [m],
    either E or l can be vectors; type get_T(\"material?\") for a list of supported materials
    """
    #get list_of supported materials from data file directory:
    #xdatafiles = [ f for f in listdir(datapath) if isfile(join(datapath,f)) ]
    xdatafiles = [ f for f in os.listdir(datapath) if os.path.isfile(os.path.join(datapath,f)) ] #OC150715
    
    name=[]
    for i in range(0, np.size(xdatafiles)):
        m=re.search('(?<=n_)\w+', xdatafiles[i])
        if m is not None:
            name.append(m.group(0))             
    
    E=np.array(E)
    l=np.array(l)
    if E.size==1 or l.size==1:
        if material in name:
            loadn=datapath+'n_'+material+'.dat'
            n=pl.loadtxt(loadn,comments='%')
            if np.min(E)>=np.min(n[:,0]) and np.max(E)<=np.max(n[:,0]):
                b=np.interp(E,n[:,0],n[:,2])
                mu=4*np.pi/get_Lambda(E,'SI')*b;
                return np.exp(-mu*l);
            else: raise xfuncs_Exception ('error: energy '+"%3.4f" %E +'[keV] out of range ('+"%3.4f" % np.min(n[:,0])+'=<E<='+"%3.4f" % np.max(n[:,0])+'eV)')
        elif material=='material?':
            #print 'list of supported materials (based on data files in directory '+datapath+':'
            #print name
            print('list of supported materials (based on data files in directory '+datapath+':')
            print(name)
        else: raise xfuncs_Exception ('error: non recognized material, please create index of refraction file first. Type "get_T?" for instructions; type get_T("material?") for list of supported materials')
    else: raise xfuncs_Exception ('error: either energy or length must be a scalar, cannot scan both energy and length at the same time.')
           
####### Help functions
#class xfuncs_Exception(Exception):
class uti_xray_mat_Exception(Exception): #OC150715
    pass
    """
    by LW 03/19/2015
    class to raise xfuncs specific exceptions
    """
