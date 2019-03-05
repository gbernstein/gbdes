#!/usr/bin/env python
"""
Process multi-extension SExtractor catalog to accomplish the following:
* Calculate median aperture corrections for stars in each chip and entire exposure,
  and put these into extension and primary headers, respectively.
* Create a header keyword MJD-MDPT if not present giving time at center of
  exposure, if one was not present.
"""
from __future__ import division,print_function
import sys
import os
import shutil
import astropy.io.fits as pf
import re
import argparse
import numpy as np
from gbutil import clippedMean
from astropy.time import Time

parser = argparse.ArgumentParser(description='Compress and process SExtractor catalogs for high-quality stars')
parser.add_argument("cat", help='Filename of catalog to augment', type=str)
parser.add_argument('--min_stars',help='Minimum number of stars to get aperture correction',
                    default=20,type=int)
parser.add_argument('--max_err',help='Max MAGERR_AUTO to use for aperture correction estimate',
                    default=0.01,type=float)
parser.add_argument('--class_star',help='Max |CLASS_STAR| to use for aperture correction estimate',
                    default=2.,type=float)

parser.add_argument('--reference_aperture',help='Which MAG_APER element is baseline',
                    default=7,type=int)
parser.add_argument('--shutter_transit_time', help='Seconds it takes for shutter to open', default=1.06, type=float)

args = parser.parse_args()

# Read input file
fits = pf.open(args.cat, 'update')
pdata = fits[0].data
phdr = fits[0].header

# These lists will hold the exposure-wide arrays of differences
# between reference aperture mag and other mags
naper = None
apcorrs = []  # For other apertures
autocorrs = [] # For MAG_AUTO
fluxrads = [] # Record FLUX_RADIUS values

# And these are chip-by-chip medians
ccd_ap09 = []
ccd_fluxrad = []

mjd = None  # MJD at start of this exposure
exptime = None  # exposure time

if "EXPTIME" in phdr:
    exptime = phdr["EXPTIME"]
if "MJD-OBS" in phdr:
    mjd = phdr["MJD-OBS"]

hduStep = 1
    
for i in range(hduStep,len(fits),hduStep):
    # Extract desirable stars from each catalog extension into new table
    hdr = fits[i].header
    data = fits[i].data

    # Get MJD if we have it
    if 'MJD-OBS' in hdr:
        if mjd is None:
            mjd = float(hdr['MJD-OBS'])
        elif mjd != float(hdr['MJD-OBS']):
            print('WARNING: disagreeing MJDs: ',mjd,hdr['MJD-OBS'])

    # and exposure time
    if 'EXPTIME' in hdr:
        if exptime is None:
            exptime = float(hdr['EXPTIME'])
        elif exptime != float(hdr['EXPTIME']):
            print('WARNING: disagreeing EXPTIMEs: ',exptime,hdr['EXPTIME'])

    # Isolate the stars we want to hold onto
    use = data['MAGERR_AUTO'] <= args.max_err
    if 'CLASS_STAR' in data.dtype.names: 
        use = np.logical_and(use, np.abs(data['CLASS_STAR'])<=args.class_star)

    tab = data[use]
    if len(tab)==0:
        continue   # Move to next header if no data here.
    
    # Calculate median aperture corrections per CCD
    # First find number of apertures used - must be the same
    # for all CCDs
    if naper is None:
        naper = len(tab['MAG_APER'][0])
        if naper <= args.reference_aperture:
            print("Aperture mag array length",naper, \
              "is too small for reference_aperture=",args.reference_aperture)
            sys.exit(1)
    else:
        if naper!=len(tab['MAG_APER'][0]):
            print("Mismatched number of photometric apertures at extension",j)
            sys.exit(1)
    # collect differences between aperture mags
    apcorr = np.array(tab['MAG_APER'])
    refmag = apcorr[:,args.reference_aperture]
    apcorr = apcorr - refmag[:,np.newaxis]
    apcorrs.append(apcorr)
    autocorrs.append(tab['MAG_AUTO'] - refmag)
    if 'FLUX_RADIUS' in tab.dtype.names:
        fluxrads.append(tab['FLUX_RADIUS'])

    if len(apcorr) >= args.min_stars:
        # Record median aperture corrections for this CCD
        ac = np.median(apcorr, axis=0)
        for i in range(naper):
            hdr['APCOR{:02d}'.format(i)] = ac[i]
        if naper>9:
            ccd_ap09.append(ac[9])
        hdr['APCORAUT'] = np.median(autocorrs[-1])
        # And median FLUX_RADIUS:
        if len(fluxrads) > 0:
            hdr['FLUXRAD'] = np.median(fluxrads[-1])
            ccd_fluxrad.append(np.median(fluxrads[-1]))

# Save exposure-wide aperture corrections
ac = np.median(np.concatenate(apcorrs), axis=0)
for i in range(naper):
    phdr['APCOR{:02d}'.format(i)] = ac[i]
phdr['APCORAUT'] = np.median(np.concatenate(autocorrs))
if len(fluxrads)>0:
    phdr['FLUXRAD'] = np.median(np.concatenate(fluxrads))
if len(ccd_fluxrad) > 20:
    mean,var,n = clippedMean(np.array(ccd_fluxrad),3.)
    phdr['FRAD_RMS'] = np.sqrt(var)
else:
    phdr['FRAD_RMS'] = 0.
if len(ccd_ap09) > 20:
    mean,var,n = clippedMean(np.array(ccd_ap09),3.)
    phdr['AP09_RMS'] = np.sqrt(var)
else:
    phdr['AP09_RMS'] = 0.

# Save the MJD at the temporal center of the exposure
if exptime is None:
    print('ERROR: add_mjdmid set but EXPTIME not found')
    sys.exit(1)
else:
    phdr['MJD-MDPT'] = mjd + 0.5*(args.shutter_transit_time + exptime)/(24.*3600.)

# Dump some basic information
print("{:6d} {:12.6f} {:5.1f} {:s} {:s} {:1s} {:+6.4f} {:6.4f} {:5.3f} {:5.3f}".format(phdr['EXPNUM'],phdr['MJD-OBS'],phdr['EXPTIME'],phdr['TELDEC'],phdr['HA'],phdr['BAND'].strip(),phdr['APCOR09'],phdr['AP09_RMS'],phdr['FLUXRAD']*2*0.264,phdr['FRAD_RMS']*2.*0.264))

fits.close()
