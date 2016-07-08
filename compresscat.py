#!/usr/bin/env python
"""
Process multi-extension SExtractor catalog to accomplish the following:
* Bring header from "LDAC_HEADER" extension into actual header of catalog extension.
* Incorporate new WCS keywords from a .head file
* Filter catalog down to clean stars with magnitude errors below selected threshold
* Calculate median aperture corrections for stars in each chip and entire exposure,
  and put these into extension and primary headers, respectively.
"""

import sys
import os
import shutil
import fitsio
import re
import argparse
import numpy as np

parser = argparse.ArgumentParser(description='Compress and process SExtractor catalogs for high-quality stars')
parser.add_argument("incat", help='Filename of input catalog', type=str)
parser.add_argument("outcat", help='Filename of output catalog', type=str)
parser.add_argument('--head', help='Name of .head file with new WCS', type=str)
parser.add_argument('--max_err',help='Maximum MAGERR_AUTO to retain',
                    default=0.01,type=float)
parser.add_argument('--imaflags_mask',help='Mask of allowed bits for IMAFLAGS_ISO',
                    default=0,type=int)
parser.add_argument('--flags_mask',help='Mask of allowed bits for FLAGS',
                    default=0,type=int)
parser.add_argument('--min_stars',help='Minimum number of stars to get aperture correction',
                    default=20,type=int)
parser.add_argument('--reference_aperture',help='Which MAG_APER element is baseline',
                    default=7,type=int)

args = parser.parse_args()

# Read input file
fitsin = fitsio.FITS(args.incat, 'r')
pdata = fitsin[0].read()
phdr = fitsin[0].read_header()

hdrs = []
tabs = []

# These lists will hold the exposure-wide arrays of differences
# between reference aperture mag and other mags
naper = None
apcorrs = []  # For other apertures
autocorrs = [] # For MAG_AUTO
psfcorrs = [] # For MAG_PSF

# These header keywords are not to be propagated from the LDAC_HEADER
skipLDAC = {'SIMPLE','BITPIX','NAXIS','NAXIS1','NAXIS2','EXTEND'}

# These header keywords are not to be propagated from the .head file
skipHead = {'FGROUPNO','ASTINST','FLXSCALE','MAGZEROP',
            'PHOTIRMS','PHOTRRMS','PHOTINST','PHOTLINK'}
    
# Open the file with new header info
if args.head:
    heads = open(args.head)
    
for i in range(2,len(fitsin),2):
    # Extract desirable stars from each catalog extension into new table
    hdr = fitsin[i].read_header()
    data = fitsin[i].read()

    # Bring LDAC headers into real headers
    ascii = fitsin[i-1][0][0][0]
    for line in ascii:
        kw = line[:8].strip()
        if len(kw)==0 or kw in skipLDAC:
            # Do not include empty or unwanted keywords
            continue
        hdr.add_record(line)
    
    if args.head:
        # Read .head files and insert into headers
        for line in heads:
            kw = line[:8].strip()
            if kw=='END':
                # We've reached the end of the header
                break
            if len(kw)==0 or kw in skipHead:
                # Do not include empty or unwanted keywords
                continue
            hdr.add_record(line.strip())

    # Isolate the stars we want to hold onto
    use = np.logical_and(np.abs(data['SPREAD_MODEL'])<=0.003,
                         data['MAGERR_AUTO'] <= args.max_err)
    use = np.logical_and(use, np.logical_not(data['FLAGS'] & args.flags_mask))
    use = np.logical_and(use, np.logical_not(data['IMAFLAGS_ISO'] & args.imaflags_mask))

    tab = data[use]

    # Calculate median aperture corrections per CCD
    # First find number of apertures used - must be the same
    # for all CCDs
    if naper is None:
        naper = len(tab['MAG_APER'][0])
        if naper <= args.reference_aperture:
            print "Aperture mag array length",naper, \
              "is too small for reference_aperture=",args.reference_aperture
            sys.exit(1)
    else:
        if naper!=len(tab['MAG_APER'][0]):
            print "Mismatched number of photometric apertures at extension",j
            sys.exit(1)
    # collect differences between aperture mags
    apcorr = tab['MAG_APER']
    refmag = apcorr[:,args.reference_aperture]
    apcorr = apcorr - refmag[:,np.newaxis]
    apcorrs.append(apcorr)
    autocorrs.append(tab['MAG_AUTO'] - refmag)
    if 'MAG_PSF' in tab.dtype.names:
        psfcorrs.append(tab['MAG_PSF'] - refmag)

    if len(apcorr) >= args.min_stars:
        # Record median aperture corrections for this CCD
        ac = np.median(apcorr, axis=0)
        for i in range(naper):
            hdr['APCOR{:02d}'.format(i)] = ac[i]
        hdr['APCORAUT'] = np.median(autocorrs[-1])
        if len(psfcorrs)>0:
            hdr['APCORPSF'] = np.median(psfcorrs[-1])

    hdrs.append(hdr)
    tabs.append(tab)

# Save exposure-wide aperture corrections
ac = np.median(np.concatenate(apcorrs), axis=0)
for i in range(naper):
    phdr['APCOR{:02d}'.format(i)] = ac[i]
    phdr['APCORAUT'] = np.median(np.concatenate(autocorrs))
    if len(psfcorrs)>0:
        phdr['APCORPSF'] = np.median(np.concatenate(psfcorrs))

# Make an output FITS and save
out = fitsio.FITS(args.outcat, 'rw', clobber=True)
for d,h in zip(tabs,hdrs):
    if 'DETPOS' in h.keys():
        out.write(d, header=h, extname=h['DETPOS'])
    else:
        out.write(d, header=h)
# Write out the primary header
out[0].write_keys(phdr)
out.close()
