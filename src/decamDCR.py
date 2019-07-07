#!/usr/bin/env python
'''
Read FOF file and generate YAML-format specifications for differential chromatic refraction terms
for each exposure (DCR).  Also calculates an improved airmass value, placing this and parallactic
angle into columns of the EXPOSURE table in the FOF file.
'''
from __future__ import division,print_function
import numpy as np
import yaml
import re
from gbutil import dmsToDegrees
import sys
import astropy.io.fits as pf
import argparse
from gbutil.decam import parallactic,dcr

# DCR amplitude in mas/mag/tan(z)
dcrConstant = {'g':45.0, 'r':8.4, 'i':3.2, 'z':1.4, 'Y':1.1}  
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
    parser.add_argument("-o", "--observatory", help='File with observatory positions', \
                        type=str)
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
    
    if args.observatory:
        observatoryTable = pf.getdata(args.observatory, 1)
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
    
    expnum_re = re.compile(r'^D(\d+)')

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
            print("Disagreeing HA's for exposure",name,ha)
        ha = dmsToDegrees(ha[0]) * 15.
        # get declination, derive parallactic angle
        dec = exposures['dec'][i]
        airmass, p = parallactic(dec, ha, lat = latitude)
        # replace airmass in table
        exposures['airmass'][i] = airmass

        # Add observatory position, if we've been given the info
        if observatoryTable is not None:
            # Find this expnum in the exposure catalog
            expnum = int(expnum_re.match(name).group(1))
            row = np.where(observatoryTable['EXPNUM']==expnum)[0]
            observatory = np.zeros(3,dtype=float)
            if len(row)<1:
                print("WARNING: exposure number",expnum,
                      "has no entry in exposure table.",row)
            elif len(row)>1:
                print("WARNING: exposure number",expnum,
                    "has multiple entries in exposure table.",row)
            else:
                observatory = observatoryTable['observatory'][row]
                #Read the observatory info and add to primary header
                while len(observatory) == 1:
                    observatory = observatory[0]
                if len(observatory)!=3:
                    print("ERROR: observatory entry in exposure table is not 3 elements")
                    sys.exit(1)
            ff['exposures'].data['OBSX'][i] = observatory[0]
            ff['exposures'].data['OBSY'][i] = observatory[1]
            ff['exposures'].data['OBSZ'][i] = observatory[2]

        # Write DCR to pixmaps
        b = bands[iinst]
        pixmaps[name+'/dcr'] = dcr(parallactic=p,
                                   airmass=airmass, band=b,
                                   referenceColor=referenceColor)

        # Write extinction gradient to pixmaps.  *Add* mags to points *closer* to zenith
        if b in nominalExtinction:
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
            print("WARNING: no nominal extinction in band",b)
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

