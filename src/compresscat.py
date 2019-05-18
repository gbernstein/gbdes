#!/usr/bin/env python
"""
Process multi-extension SExtractor catalog to accomplish the following:
* Bring header from "LDAC_HEADER" extension into actual header of catalog extension.
* Incorporate new WCS keywords from a .head file
* Filter catalog down to clean stars with magnitude errors below selected threshold
* Calculate median aperture corrections for stars in each chip and entire exposure,
  and put these into extension and primary headers, respectively.
* Choose a star flat calibration epoch to apply to this exposure (CALEPOCH) based on its MJD.
* Create a header keyword MJD-MDPT if not present giving time at center of
  exposure, if one was not present.
"""
from __future__ import division,print_function
import sys
import os
import shutil
import fitsio
import re
import argparse
import numpy as np
from gbutil import clippedMean
from astropy.time import Time
from gbutil.decam import EpochFinder

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
parser.add_argument('--add_mjdmid', help='Add MJD-MDPT to primary header', default=True, type=bool)
parser.add_argument('--shutter_transit_time', help='Seconds it takes for shutter to open', default=1.06, type=float)

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
fluxrads = [] # Record FLUX_RADIUS values

# And these are chip-by-chip medians
ccd_ap09 = []
ccd_fluxrad = []

# These header keywords are not to be propagated from the LDAC_HEADER
skipLDAC = {'SIMPLE','BITPIX','NAXIS','NAXIS1','NAXIS2','EXTEND','CONTINUE'}

# These header keywords are not to be propagated from the .head file
skipHead = {'FGROUPNO','ASTINST','FLXSCALE','MAGZEROP',
            'PHOTIRMS','PHOTRRMS','PHOTINST','PHOTLINK'}
    
mjd = None  # MJD at start of this exposure
exptime = None  # exposure time

ef = EpochFinder()

# Open the file with new header info
if args.head:
    heads = open(args.head)
    
# Is this an LDAC-format file, or are headers already in the data?
h1 = fitsin[1].read_header()
isLDAC = 'EXTNAME' in h1 and \
         (h1['EXTNAME']=='LDAC_HEADER' or h1['EXTNAME']=='LDAC_IMHEAD')

if isLDAC:
    hduStep = 2
else:
    hduStep = 1
    
for i in range(hduStep,len(fitsin),hduStep):
    # Extract desirable stars from each catalog extension into new table
    hdr = fitsin[i].read_header()
    data = fitsin[i].read()

    if isLDAC:
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
    if 'SPREAD_MODEL' in data:  ## ??? string comparison to numpy array potential problem ???
        use = np.logical_and(use, np.abs(data['SPREAD_MODEL'])<=0.003)
    use = np.logical_and(use, np.logical_not(data['FLAGS'] & args.flags_mask))
    use = np.logical_and(use, np.logical_not(data['IMAFLAGS_ISO'] & args.imaflags_mask))

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
    apcorr = tab['MAG_APER']
    refmag = apcorr[:,args.reference_aperture]
    apcorr = apcorr - refmag[:,np.newaxis]
    apcorrs.append(apcorr)
    autocorrs.append(tab['MAG_AUTO'] - refmag)
    if 'MAG_PSF' in tab.dtype.names:
        psfcorrs.append(tab['MAG_PSF'] - refmag)
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
        if len(psfcorrs)>0:
            hdr['APCORPSF'] = np.median(psfcorrs[-1])
        # And median FLUX_RADIUS:
        if len(fluxrads) > 0:
            hdr['FLUXRAD'] = np.median(fluxrads[-1])
            ccd_fluxrad.append(np.median(fluxrads[-1]))
    hdrs.append(hdr)
    tabs.append(tab)

# Save exposure-wide aperture corrections
ac = np.median(np.concatenate(apcorrs), axis=0)
for i in range(naper):
    phdr['APCOR{:02d}'.format(i)] = ac[i]
phdr['APCORAUT'] = np.median(np.concatenate(autocorrs))
if len(psfcorrs)>0:
    phdr['APCORPSF'] = np.median(np.concatenate(psfcorrs))
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

# Save calibration epoch
epoch =  ef(mjd)
phdr['CALEPOCH'] =epoch    
# Save the MJD at the temporal center of the exposure
if args.add_mjdmid:
 if exptime is None:
     print('ERROR: add_mjdmid set but EXPTIME not found')
     sys.exit(1)
 else:
    phdr['MJD-MDPT'] = mjd + 0.5*(args.shutter_transit_time + exptime)/(24.*3600.)

# Dump some basic information
print("{:6d} {:12.6f} {:5.1f} {:s} {:s} {:1s} {:+6.4f} {:6.4f} {:5.3f} {:5.3f}".format(hdr['EXPNUM'],hdr['MJD-OBS'],hdr['EXPTIME'],hdr['TELDEC'],hdr['HA'],hdr['BAND'].strip(),phdr['APCOR09'],phdr['AP09_RMS'],phdr['FLUXRAD']*2*0.264,phdr['FRAD_RMS']*2.*0.264), epoch)

# Make an output FITS and save
out = fitsio.FITS(args.outcat, 'rw', clobber=True)
for d,h in zip(tabs,hdrs):
    if 'DETPOS' in h:
        out.write(d, header=h, extname=h['DETPOS'])
    else:
        out.write(d, header=h)
# Write out the primary header
out[0].write_keys(phdr)
out.close()
