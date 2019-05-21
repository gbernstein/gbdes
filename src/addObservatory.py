#!/usr/bin/env python
from __future__ import division,print_function

import sys
import argparse
import astropy.io.fits as pf
import numpy as np

parser = argparse.ArgumentParser(description=\
            'Add observatory position and mjd-mid to headers of FITS catalog files')
parser.add_argument("exposureTable", help='filename of FITS table of exposure info', type=str)
parser.add_argument("catalog", help='Catalogs to update',type=str,nargs="+")

def addObservatory(exposureTable, fitsFile):
    ''' Find first extension having an EXPNUM keyword. Use it to find observatory posn
    and add this to the primary header.
    '''

    with pf.open(fitsFile,'update') as ff:
        expnum = None
        for hdu in ff:
            if 'EXPNUM' in hdu.header:
                expnum = hdu.header['EXPNUM']
                break
        if expnum is None:
            print("WARNING: no EXPNUM found for catalog",fitsFile)
            return
        # Find this expnum in the exposure catalog
        row = np.where(exposureTable['expnum']==expnum)[0]
        if len(row)<1:
            print("WARNING: exposure number",expnum,"from file",fitsFile,
                         "has no entry in exposure table.")
            return
        elif len(row)>1:
            print("WARNING: exposure number",expnum,"from file",fitsFile,
                    "has multiple entries in exposure table.")
            return

        #Read the observatory info and add to primary header
        observatory = exposureTable['observatory'][row]
        while len(observatory) == 1:
            observatory = observatory[0]
        mjd = exposureTable[row]['mjd_mid'][0]
        if len(observatory)!=3:
            print("ERROR: observatory entry in exposure table is not 3 elements")
            sys.exit(1)
        ff[0].header['OBSX'] = observatory[0]
        ff[0].header['OBSY'] = observatory[1]
        ff[0].header['OBSZ'] = observatory[2]
        ff[0].header['MJD-MID'] = mjd
if __name__=="__main__":
    args = parser.parse_args()
    exposures = pf.getdata(args.exposureTable,1)
    for n in args.catalog:
        print(n)
        addObservatory(exposures, n)
    sys.exit(0)
