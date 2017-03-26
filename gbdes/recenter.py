#!/usr/bin/env python
"""
Program to replace the nominal RA, Dec of DECam exposures with the average of the centers of the
S4 and N4 CCDs in the config or fof files.  The WCSIN headers are to be interpreted by the astropy
wcs routines.
"""

import sys
import numpy as np
import astropy.io.fits as pf
import astropy.wcs as wcs
from gmbpy import headerFromString

def fixRADec(fits):
    try:
        exposures = fits['Exposures'].data
        extensions = fits['Extensions'].data
    except KeyError:
        print "File does not have necessary extensions"
        return 1
    for iexpo,name in enumerate(exposures['name']):
        print "Exposure",name
        iinst = exposures['instrumentnumber'][iexpo]
        if iinst<0:
            # Not using a decam catalog, so skip
            continue
        # Find the instrument with this number and get
        # device numbers for N4, S4
        in4 = None
        is4 = None
        for h in fits[1:]:
            if h.header['extname'].lower()!='Instrument'.lower() or h.header['number']!=iinst:
                continue
            wn = np.where(h.data['name']=='N4')[0]
            ws = np.where(h.data['name']=='S4')[0]
            if len(wn)==1 and len(ws)==1:
                in4 = wn[0]
                is4 = ws[0]
        if in4 is None or is4 is None:
            print "Warning: Did not find N4 and S4 devices for instrument number",iinst
            continue

        en4 = np.logical_and(extensions['exposure']==iexpo, extensions['device']==in4)
        es4 = np.logical_and(extensions['exposure']==iexpo, extensions['device']==is4)
        if np.count_nonzero(en4)!=1 or np.count_nonzero(es4)!=1:
            print "Warning: did not find N4 and S4 extensions for exposure",name
            continue
        w = wcs.WCS(headerFromString(extensions['WCSIN'][en4][0].split('\n')))
        radecn = w.wcs_pix2world( [[1024,2048]],1)
        w = wcs.WCS(headerFromString(extensions['WCSIN'][es4][0].split('\n')))
        radecs = w.wcs_pix2world( [[1024,2048]],1)
        radec = 0.5*(radecn+radecs)  # Could check for wrap around of RA???
        ra = radec.flatten()[0]
        dec = radec.flatten()[1]
        print " old:",exposures['RA'][iexpo],exposures['Dec'][iexpo]
        print " new:",ra,dec
        exposures['RA'][iexpo] = ra
        exposures['Dec'][iexpo] = dec
    return 0

def fieldCoords(fits):
    """
    Update or create a new column in the exposure table giving the exposure
    axis position in gnomonic position about its field center.
    """
    # check for field and instrument exposures
    try:
        exposures = fits['Exposures'].data
        fields = fits['Fields'].data
    except KeyError:
        print "File does not have necessary extensions"
        return 1
    # Add columns to the exposure table ???
    for iexpo,name in enumerate(exposures['name']):
        print "Exposure",name
        iinst = exposures['instrumentnumber'][iexpo]
        if iinst<0:
            # Not using a decam catalog, so skip
            continue
        #  find field RA & dec
        ifield = exposures['FIELDNUMBER'][iexpo]
        fra = fields['RA'][ifield]
        fdec = fields['Dec'][ifield]
        ra = exposures['RA'][iexpo]
        dec = exposures['Dec'][iexpo]
        # Calculate gnomonic projection (degrees here by default)
        w = wcs.WCS(naxis=2)
        w.wcs.ctype = ["RA---TAN","DEC--TAN"]
        w.wcs.crval = np.array([fra, fdec])
        dx, dy = w.wcs_world2pix(ra,dec,1)
        print ra, fra, dx
        print dec,fdec, dy
        # save dx[0], dy[0]
    return 0

if __name__=='__main__':
    if len(sys.argv) != 2:
        print "Usage: recenter.py filename"
        print " <filename> is the output of configure.py that will have exposure"
        print "            RA/Dec values updated."
        sys.exit(1)
    filename = sys.argv[1]
    fits = pf.open(filename, 'update')
    code = fixRADec(fits)
    ##if code==0:
    ##    code = fieldCoords(fits)
    fits.close()
    sys.exit(code)

        
        
