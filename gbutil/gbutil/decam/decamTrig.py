#!/usr/bin/env python
'''
Routines for some coordinate calculations used in DECam and elsewhere
'''
import numpy as np
import sys
import astropy.io.fits as pf
import astropy.wcs as wcs

def parallactic(dec, ha, lat=-30.1716):
    '''Function will calculate airmass and
    parallactic angle (in radians) given
    input of source declination, hour angle,
    and observatory latitude (in degrees)
    '''
    # First get cosine of zenith angle
    dtor = np.pi / 180.
    d = dec*dtor
    if (ha>180.):
        h = (ha-360.)*dtor  # Make negative HA instead of large
    else:
        h = ha*dtor
    l = lat*dtor
    cosz = np.sin(l)*np.sin(d) + np.cos(l)*np.cos(d)*np.cos(h)
    sinz = np.sqrt(1-cosz*cosz)  # zenith angle always 0-180
    airmass = 1./cosz
    # Now the parallactic angle
    if sinz<=0:
        # at (anti-)zenith already, p is undefined, return 0
        return airmass,0.
    cosp = (np.sin(l)*np.cos(d)-np.cos(l)*np.sin(d)*np.cos(h))/sinz
    if np.abs(cosp)>1.:
        cosp = np.sign(cosp)
    p = np.arccos(cosp)* np.sign(h)
    return airmass,p

def dcr(parallactic, airmass, band, referenceColor=0.61):
    '''Return a dict specifying a valid PixelMap for the 
    differential chromatic refraction in this exposure
    with DECam.
    The referenceColor is the g-i for which DCR is defined as zero.
    '''
    # DCR amplitude in mas/mag/tan(z)
    dcrConstant = {'g':45.0, 'r':8.4, 'i':3.2, 'z':1.4, 'Y':1.1}  
    if band in dcrConstant:
        # Convert DCR amplitude from mas to degrees
        ampl = dcrConstant[band] / (3600.*1000.)
        dcry = ampl * np.cos(parallactic)*np.sqrt(airmass*airmass-1)
        dcrx = ampl * np.sin(parallactic)*np.sqrt(airmass*airmass-1)
        # Specify the distortion
        pixmap = {'Type':'Color',
                  'Reference':referenceColor,
                  'Function':{'Type':'Constant',
                              'Parameters':[float(dcrx),float(dcry)]}}
    else:
        # not a band that needs DCR
        pixmap = {'Type':'Identity'}

    return pixmap

def recenter(hn4, hs4):
    '''Return new optic-axis RA, Dec (in degrees) given the FITS WCS
    headers for CCDs N4 and S4.  Define axis as midpoint of these.
    '''
    w = wcs.WCS(hn4)
    radecn = w.wcs_pix2world( [[1024,2048]],1)
    w = wcs.WCS(hs4)
    radecs = w.wcs_pix2world( [[1024,2048]],1)
    radec = 0.5*(radecn+radecs)  # Could check for wrap around of RA???
    ra = radec.flatten()[0]
    dec = radec.flatten()[1]
    return ra, dec


