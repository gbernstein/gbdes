#!/usr/bin/env python
'''
Read FOF file and generate YAML-format specifications for differential chromatic refraction terms
for each exposure (DCR).  Also calculates an improved airmass value, placing this and parallactic
angle into columns of the EXPOSURE table in the FOF file.
'''
import numpy as np
import yaml
from gbutil import dmsToDegrees
import sys
import astropy.io.fits as pf
import argparse



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

dcrConstant = {'g':45.0, 'r':8.4}  # DCR amplitude in mas/mag/tan(z)
latitude = -30.1716  # Observatory latitude, degrees
referenceColor = 0.61  # zeropoint for color terms
# Extinction coeffs from https://cdcvs.fnal.gov/redmine/projects/descalibration/wiki
nominalExtinction = {'u':0.45, 'g':0.20,'r':0.10,'i':0.07,'z':0.08,'Y':0.07}

if __name__=='__main__':
    parser = argparse.ArgumentParser(description=\
     'Create DCR pixelmaps for exposures in FOF file and update airmass entries')
    parser.add_argument("fof", help='Filename of FOF file to update', type=str)
    parser.add_argument("dcr", help='Filename of output YAML DCR specs', type=str)
    parser.add_argument("gradient", help='Filename of output YAML extinction gradients', type=str)
    parser.add_argument('--no_update',help='Do not update FOF file airmasses',
                        action='store_true')

    args = parser.parse_args()
    fits = args.fof
    dcrFile = args.dcr
    gradientFile = args.gradient

    pixmaps = {}
    photomaps = {}
    pixmaps['PixelMapCollection'] = ' '.join(sys.argv)
    photomaps['PhotoMapCollection'] = ' '.join(sys.argv)
    
    if args.no_update:
        ff = pf.open(fits)
    else:
        ff = pf.open(fits,'update')

    # Get bands of all instruments, assign band to exposures
    bands = {}
    for hdu in ff[1:]:
        if hdu.header['EXTNAME'].lower()=='Instrument'.lower():
            bands[hdu.header['NUMBER']] = hdu.header['BAND']

    exposures = ff['exposures'].data
    extns = ff['extensions'].data
    
    # Read exposure table
    for i,expo in enumerate(exposures['name']):
        iinst = exposures['instrumentnumber'][i]
        if iinst < 0:
            # Not a real instrument, no changes here
            continue
        name = "".join(expo)
        # Find all extensions with this exposure
        used = extns['exposure']==i
        # Get HA from them, do they agree?
        ha = np.unique(extns['HA'][used])
        if len(ha) > 1:
            print "Disagreeing HA's for exposure",name,ha
        ha = dmsToDegrees(ha[0]) * 15.
        # get declination, derive parallactic angle
        dec = exposures['dec'][i]
        airmass, p = parallactic(dec, ha, lat = latitude)
        # replace airmass in table
        exposures['airmass'][i] = airmass

        # Write DCR to pixmaps
        b = bands[iinst]
        if b in dcrConstant.keys():
            # Convert DCR amplitude from mas to degrees
            ampl = dcrConstant[b] / (3600.*1000.)
            dcry = ampl * np.cos(p)*np.sqrt(airmass*airmass-1)
            dcrx = ampl * np.sin(p)*np.sqrt(airmass*airmass-1)
            # Specify the distortion
            pixmap = {'Type':'Color',
                      'Reference':referenceColor,
                      'Function':{'Type':'Constant',
                                  'Parameters':[float(dcrx),float(dcry)]}}
            pixmaps[name+'/dcr'] = pixmap
        else:
            # not a band that needs DCR
            pixmaps[name+'/dcr'] = {'Type':'Identity'}

        # Write extinction gradient to pixmaps.  *Add* mags to points *closer* to zenith
        if b in nominalExtinction.keys():
            # Get local extinction gradient in mag per degree
            ampl = nominalExtinction[b] * airmass * np.sqrt(airmass**2-1) * np.pi/180.;
            dcry = ampl * np.cos(p) 
            dcrx = ampl * np.sin(p)
            # Specify the distortion.  XMin/Max will default to +-1 to give native scaling
            photomap = {'Type':'Poly',
                        'UsePixelCoords':False,
                        'Poly':{'OrderX':1,
                                'SumOrder': True,
                                'Coefficients':[0.,float(dcrx),float(dcry)]} }
            photomaps[name+'/gradient'] = photomap
        else:
            # Extinction not known
            print "WARNING: no nominal extinction in band",b
            photomaps[name+'/gradient'] = {'Type':'Identity'}
            
    # write yaml
    fout = open(dcrFile,'w')
    yaml.dump(pixmaps,fout)
    fout.close()
    fout = open(gradientFile,'w')
    yaml.dump(photomaps,fout)
    fout.close()
    ff.close()
    sys.exit(0)

