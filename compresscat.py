#!/usr/bin/env python
"""
Process multi-extension SExtractor catalog to accomplish the following:
* Bring header from "LDAC_HEADER" extension into actual header of catalog extension.
* Incorporate new WCS keywords from a .head file
* Filter catalog down to clean stars with magnitude errors below selected threshold
* Calculate median aperture corrections for stars in each chip and entire exposure,
  and put these into extension and primary headers, respectively.
* Choose a star flat calibration epoch to apply to this exposure (CALEPOCH) based on its MJD.
"""

import sys
import os
import shutil
import fitsio
import re
import argparse
import numpy as np
from gmbpy import clippedMean
from astropy.time import Time

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
fluxrads = [] # Record FLUX_RADIUS values

# And these are chip-by-chip medians
ccd_ap09 = []
ccd_fluxrad = []

# These header keywords are not to be propagated from the LDAC_HEADER
skipLDAC = {'SIMPLE','BITPIX','NAXIS','NAXIS1','NAXIS2','EXTEND','CONTINUE'}

# These header keywords are not to be propagated from the .head file
skipHead = {'FGROUPNO','ASTINST','FLXSCALE','MAGZEROP',
            'PHOTIRMS','PHOTRRMS','PHOTINST','PHOTLINK'}
    
mjd = None  # MJD of this exposure
def mjdOfEpoch(mjd):
    # Return mjd of date specified by 8-character epoch
    return Time(epoch[:4]+'-'+epoch[4:6]+'-'+epoch[6:8], 
                format='fits',scale='utc').mjd
class EpochFinder(object):
    '''Function class which returns the star flat epoch nearest to specifed input MJD
    that does not have an intervening camera event.  Returns '00000000' if no star flats
    occur in the same interval between events.
    '''
    # Epochs of star flats and of camera "events" when calibration changes.
    sfEpochs = ['20121223','20130221','20130829','20131115','20140118','20140807','20141105',
                '20150204','20150926', '20160209', '20160223','20160816','20161117']
    warmups = ['20121230','20130512','20130722','20131015','20140512','20141201','20150625',
                '20150725','20150809','20150825','20160219','20161013']
    cooldowns=['20121226','20151126']
    nogood = '00000000'
    def __init__(self):
        self.sfMjds = np.array([mjdOfEpoch(e) for e in self.sfEpochs])
        self.eventMjds = np.array([mjdOfEpoch(e) for e in self.warmups + self.cooldowns].sort())

        return
    def __call__(self, mjd):
        if mjd is None:
            return nogood
        # which events are before, after our mjd?
        before = mjd >= self.eventMjds
        # Mark which star flat MJDs are in same interval between events
        if not np.any(before):
            # Our mjd is before any events
            same = sfMjds < eventMjds[0]
        elif np.all(before):
            # Our mjd is after all events
            same = sfMjds >= eventMjds[-1]
        else:
            # Our mjd is between two events, get index of preceding one
            precede = np.where(before)[0][-1]
            same = np.logical_and(sfMjds>=eventMjds[precede],
                                  sfMjds < eventMjds[precede+1])
        
        if not same.any():
            # No star flats in the same event interval.
            return nogood
        sameIndices = np.where(same)[0]
        closest = np.argmin(np.abs(sfMjds[same]-mjd))
        return sfEpochs[sameIndices[closest]]

ef = EpochFinder()

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

    # Get MJD if we have it
    if 'MJD-OBS' in hdr.keys():
        if mjd is not None:
            mjd = float(hdr['MJD-OBS'])
        elif mjd != float(hdr['MJD-OBS']):
            print 'WARNING: disagreeing MJDs: ',mjd,hdr['MJD-OBS']

    # Isolate the stars we want to hold onto
    use = np.logical_and(np.abs(data['SPREAD_MODEL'])<=0.003,
                         data['MAGERR_AUTO'] <= args.max_err)
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

# Dump some basic information
print "{:6d} {:12.6f} {:5.1f} {:s} {:s} {:1s} {:+6.4f} {:6.4f} {:5.3f} {:5.3f}".format(hdr['EXPNUM'],hdr['MJD-OBS'],hdr['EXPTIME'],hdr['TELDEC'],hdr['HA'],hdr['BAND'].strip(),phdr['APCOR09'],phdr['AP09_RMS'],phdr['FLUXRAD']*2*0.264,phdr['FRAD_RMS']*2.*0.264), 

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
