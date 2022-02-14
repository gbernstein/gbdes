""" Generically useful routines
"""
from __future__ import division,print_function
import astropy.io.fits as pf
import numpy as np
import sys
import re

def readLDACHeader(filename, extn):
    """ Return a pyfits Header object from LDAC format"""
    ascii = pf.getdata(filename, extn)[0][0]
    return headerFromString(ascii)

def headerFromString(ascii):
    """ Return pyfits header from either an array of strings, or a long string with
    80-char "cards" as per FITS std
    """
    cardlist = []
    # Note newer PyFITS would let us pf.append(pf.Card.fromstring(i))
    if type(ascii) == str:
        for j in range(0, len(ascii), 80):
            cardlist.append(pf.Card.fromstring(ascii[j:j+80]))
    else:
        for i in ascii:
            cardlist.append(pf.Card.fromstring(i))
    return pf.Header(cardlist)

def mergeLDACs(infiles, outfile, keywordForExtname=None):
    """ Merge the LDAC format tables in infiles into a single
    output FITS files.  Put the LDAC header back into the true
    header of the output HDU.  Give new extension a name from the
    value associated with chosen keyword, if desired.
    """
    outlist = pf.HDUList([pf.PrimaryHDU()])

    for infits in infiles:
        f = pf.open(infits)
        hdunum = 1
        while hdunum < len(f):
            if f[hdunum].header['EXTNAME']=='LDAC_IMHEAD':
                # Have an LDAC here.  Use it.
                h = readLDACHeader(infits, hdunum)
                hdunum += 1
                if hdunum >= len(f) or f[hdunum].header['EXTNAME']!='LDAC_OBJECTS':
                    print('Missing/wrong LDAC_OBJECTS in ',infits,\
                          ' at extension ',hdunum)
                    sys.exit(1);
                hdu = f[hdunum]
                hdu.header.extend(h)
                if (keywordForExtname!=None):
                    hdu.header['EXTNAME'] = hdu.header[keywordForExtname]
                outlist.append(hdu)
            hdunum += 1
    outlist.writeto(outfile)
    return


def dmsToDegrees(dms):
    """ Convert colon-separted degrees/min/sec string to decimal degrees"""
    ldms = dms.split(':')
    minus = ldms[0].strip()[0]=='-'
    degrees = float(ldms[0])
    if degrees<0.: degrees *= -1
    minutes = 0.
    if len(ldms)>2: minutes = float(ldms[2])/60.
    if len(ldms)>1: minutes += float(ldms[1])
    degrees += minutes/60.
    if minus: degrees *= -1.
    return degrees

def isInt(s):
    """ Return true if input string is a valid integer representation"""
    try:
        int(s)
        return True
    except ValueError:
        return False

def isFloat(s):
    """ Return true if input string is a valid float representation"""
    try:
        float(s)
        return True
    except ValueError:
        return False

defaultReferenceDay=(2012,11,1)
""" Default start time for our day clock is 0h UT on Nov 1 2012"""

def days(ymdhms, reference=defaultReferenceDay):
    """ Return number of days since midnight of reference date, given yr/mo/day/hrs/min/sec

    Input should be tuple of integers or floats or strings that convert to such.
    Can truncate the tuple at any point after the day and assume remainder are zeros.
    The reference field gives the zeropoint for data counting.
    Simple routine is not going to work for leap years or leap seconds.  And it's slow.
    """
    daysInMonth=31,28,31,30,31,30,31,31,30,31,30,31
    initDay=[0,0]
    for i in daysInMonth:
        initDay.append(initDay[-1] + i)

    input = [0, 0, 0, 0, 0, 0.]
    ref = [0, 0, 0, 0, 0, 0.]
    ref[0:len(reference)] = reference
    input[0:len(ymdhms)] = ymdhms

    days = float(input[5]) - float(ref[5])
    days = days/60. + float(input[4]) - float(ref[4])
    days = days/60. + float(input[3]) - float(ref[3])
    days = days/24. + float(input[2]) + initDay[int(input[1])]
    days -= float(ref[2]) + initDay[int(ref[1])]
    days += (int(input[0])-int(ref[0]))*initDay[13]
    return days

def dateobsToDay(dateobs,reference=defaultReferenceDay):
    """Convert YYYY-MM-DDTHH:MM:SS.SSSS format to float days since reference
    will split at space if there is no T between date and time.
    And will get rid of the +00:00 at end if it is there. 
    """
    dt = dateobs.split('T')
    if len(dt)==1:
        # Split at white space if there was no T:
        dt = dateobs.split()
    return days( dt[0].split('-') + dt[1].split('+')[0].split(':'), reference)

def parallacticAngle(zenithDistance, hourAngle, declination, latitude=-(30+10./60.)):
    """Return parallactic angle in degrees from N through E.  Inputs in degrees."""
    dtor = 3.14159265/180.
    parallactic = np.arcsin(np.cos(latitude*dtor) * np.sin(hourAngle*dtor) / 
                            np.sin(zenithDistance*dtor)) / dtor
    if np.sin(latitude*dtor) > np.cos(zenithDistance*dtor)*np.sin(declination*dtor):
        parallactic = 180. - parallactic
    return parallactic
    
def sectionToSlices(sec):
    """ Convert a pixel range in FITS/IRAF string format '[x0:x1, y0:y1, ...]' to
    a tuple of slices (....,sy,sx) giving numpy index ranges.
    including shift from numpy 0-indexing to FITS 1-indexing, and also
    reversal of index order.
    ***Not prepared to deal with wildcards in image sections right now***
    """
    out = ()
    ranges = sec.strip('[]').split(',')
    for r in reversed(ranges):
        r2 = r.split(':')
        if not len(r2)==2:
            print('Bad range in sectionToRange: ',r)
            return
        s = slice( int(r2[0])-1, int(r2[1]) )
        out += (s,)
    return out

def slicesToSection(slices):
    """ Convert a pixel range in given by tuple of slices to a FITS/IRAF
    image section specification string '[x0:x1,y0:y1,...]'
    including shift from numpy 0-indexing to FITS 1-indexing, and also
    reversal of index order.
    """
    out = '['
    for s in reversed(slices):
        out += '{:d}:{:d},'.format(s.start+1, s.stop)
    # Now strip the trailing , and add closing bracket:
    return out[:-1]+']'

def extractWCS(headers):
    """
    Extract from a FITS header the keywords that specify the WCS.
    """
    wcsre = re.compile(r"^(CR(PIX|VAL)[12]|PV[12]_\d*|CD[12]_[12]|C(TYPE|UNIT)[12]|EQUINOX|RADESYS)\s*$")

    if type(headers)==pf.Header or type(headers)==dict:
        hlist = [headers]
    else:
        hlist = headers
    newhdr = {}
    for h in hlist:
        for k,v in h.items():
            if wcsre.match(k):
                newhdr[k] = v
    out = pf.Header()
    for k,v in newhdr.items():
        out.append( (k,v))

    return out



    
