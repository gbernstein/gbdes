""" Routines here that give or manipulate images with formats specific to DECam.
"""
from __future__ import division,print_function
import numpy as np
import astropy.io.fits as pf
import re

plotColors = { 'g':'g','r':'r','i':'c','z':'m','Y':'y'}

amps = ('A','B') # Possible amplifier choices
dps = [] # List of detpos values in "pretty" order
for i in range(1,32): dps.append('N'+str(i))
for i in range(1,32): dps.append('S'+str(i))


def overscanColumns(detpos, amp):
    """ Returns source, target where each is a slice giving columns
    that contain overscan and from which overscan should be subtracted, respectively 
    """
    left = slice(10,50),slice(56,1080)
    right= slice(2110,2150),slice(1080,2104)
    if detpos[0]=='N':
        if amp=='A': return left
        else: return right
    elif detpos[0]=='S':
        if amp=='A': return right
        else: return left
    else:
        raise Exception('Gave invalid detpos to overscanColumns: '+str(detpos))
    return

def centralAmpPix(detpos, amp):
    """ Return a 2d slice tuple specifying central 512x512 region (1/8 of pix) of an amp
    after trimming to 2k x 4k has been done
    """
    left = ( slice(1792,2304), slice(256,768) )
    right =( slice(1792,2304), slice(1280,1792) )
    if detpos[0]=='N':
        if amp=='A': return left
        else: return right
    elif detpos[0]=='S':
        if amp=='A': return right
        else: return left
    else:
        raise Exception('Gave invalid detpos to centralAmpPix: '+str(detpos))
    return

def allAmpPix(detpos, amp):
    """ Return a 2d slice tuple specifying all of the pixels that are from amp A, after trimming.
    """
    left = slice(None),slice(None,1024)
    right = slice(None),slice(1024,None)
    if detpos[0]=='N':
        if amp=='A': return left
        else: return right
    elif detpos[0]=='S':
        if amp=='A': return right
        else: return left
    else:
        raise Exception('Gave invalid detpos to allAmpPix: '+str(detpos))

        
ccdnums =  {'S29':  1,
            'S30':  2,
            'S31':  3,
            'S25':  4,
            'S26':  5,
            'S27':  6,
            'S28':  7,
            'S20':  8,
            'S21':  9,
            'S22':  10,
            'S23':  11,
            'S24':  12,
            'S14':  13,
            'S15':  14,
            'S16':  15,
            'S17':  16,
            'S18':  17,
            'S19':  18,
            'S8':  19,
            'S9':  20,
            'S10':  21,
            'S11':  22,
            'S12':  23,
            'S13':  24,
            'S1':  25,
            'S2':  26,
            'S3':  27,
            'S4':  28,
            'S5':  29,
            'S6':  30,
            'S7':  31,
            'N1':  32,
            'N2':  33,
            'N3':  34,
            'N4':  35,
            'N5':  36,
            'N6':  37,
            'N7':  38,
            'N8':  39,
            'N9':  40,
            'N10':  41,
            'N11':  42,
            'N12':  43,
            'N13':  44,
            'N14':  45,
            'N15':  46,
            'N16':  47,
            'N17':  48,
            'N18':  49,
            'N19':  50,
            'N20':  51,
            'N21':  52,
            'N22':  53,
            'N23':  54,
            'N24':  55,
            'N25':  56,
            'N26':  57,
            'N27':  58,
            'N28':  59,
            'N29':  60,
            'N30':  61,
            'N31':  62}

# DECam CCD corners in system where they are abutted into one big array
# Note that x coord is offset by 2048 because of guide CCDs.
ccdCorners = {'N1': ( 14337,     1),
              'N2': ( 14337,  4097),
              'N3': ( 14337,  8193),
              'N4': ( 14337, 12289),
              'N5': ( 14337, 16385),
              'N6': ( 14337, 20481),
              'N7': ( 14337, 24577),
              'N8': ( 16385,  2049),
              'N9': ( 16385,  6145),
              'N10': ( 16385, 10241),
              'N11': ( 16385, 14337),
              'N12': ( 16385, 18433),
              'N13': ( 16385, 22529),
              'N14': ( 18433,  2049),
              'N15': ( 18433,  6145),
              'N16': ( 18433, 10241),
              'N17': ( 18433, 14337),
              'N18': ( 18433, 18433),
              'N19': ( 18433, 22529),
              'N20': ( 20481,  4097),
              'N21': ( 20481,  8193),
              'N22': ( 20481, 12289),
              'N23': ( 20481, 16385),
              'N24': ( 20481, 20481),
              'N25': ( 22529,  6145),
              'N26': ( 22529, 10241),
              'N27': ( 22529, 14337),
              'N28': ( 22529, 18433),
              'N29': ( 24577,  8193),
              'N30': ( 24577, 12289), ###  dead CCD!!
              'N31': ( 24577, 16385),
              'S1': ( 12289,     1),
              'S2': ( 12289,  4097),
              'S3': ( 12289,  8193),
              'S4': ( 12289, 12289),
              'S5': ( 12289, 16385),
              'S6': ( 12289, 20481),
              'S7': ( 12289, 24577),
              'S8': ( 10241,  2049),
              'S9': ( 10241,  6145),
              'S10': ( 10241, 10241),
              'S11': ( 10241, 14337),
              'S12': ( 10241, 18433),
              'S13': ( 10241, 22529),
              'S14': (  8193,  2049),
              'S15': (  8193,  6145),
              'S16': (  8193, 10241),
              'S17': (  8193, 14337),
              'S18': (  8193, 18433),
              'S19': (  8193, 22529),
              'S20': (  6145,  4097),
              'S21': (  6145,  8193),
              'S22': (  6145, 12289),
              'S23': (  6145, 16385),
              'S24': (  6145, 20481),
              'S25': (  4097,  6145),
              'S26': (  4097, 10241),
              'S27': (  4097, 14337),
              'S28': (  4097, 18433),
              'S29': (  2049,  8193),
              'S30': (  2049, 12289),
              'S31': (  2049, 16385) }

""" Next is a dictionary giving the approx corner positions of each CCD on the sky, in degrees
from center of the focal plane.  Tuple is (xmin, xmax, ymin, ymax) with x to E and y to N.
"""
ccdBounds = {'N1': (-1.0811, -0.782681, -0.157306, -0.00750506),
             'N2': (-0.771362, -0.472493, -0.157385, -0.00749848), 
             'N3': (-0.461205, -0.161464, -0.157448, -0.00749265), 
             'N4': (-0.150127, 0.149894, -0.15747, -0.00749085), 
             'N5': (0.161033, 0.460796, -0.157638, -0.0074294), 
             'N6': (0.472171, 0.771045, -0.157286, -0.00740563), 
             'N7': (0.782398, 1.08083, -0.157141, -0.0074798), 
             'N8': (-0.92615, -0.627492, -0.321782, -0.172004), 
             'N9': (-0.616455, -0.317043, -0.322077, -0.172189), 
             'N10': (-0.305679, -0.00571999, -0.322071, -0.17217), 
             'N11': (0.00565427, 0.305554, -0.322243, -0.172254), 
             'N12': (0.31684, 0.616183, -0.322099, -0.172063), 
             'N13': (0.627264, 0.925858, -0.321792, -0.171887), 
             'N14': (-0.926057, -0.62726, -0.485961, -0.336213), 
             'N15': (-0.616498, -0.317089, -0.486444, -0.336606), 
             'N16': (-0.30558, -0.00578257, -0.486753, -0.336864), 
             'N17': (0.00532179, 0.305123, -0.486814, -0.33687), 
             'N18': (0.316662, 0.616018, -0.486495, -0.336537), 
             'N19': (0.62708, 0.92578, -0.485992, -0.336061), 
             'N20': (-0.770814, -0.471826, -0.650617, -0.500679), 
             'N21': (-0.460777, -0.161224, -0.650817, -0.501097), 
             'N22': (-0.149847, 0.149886, -0.650816, -0.501308), 
             'N23': (0.161001, 0.460566, -0.650946, -0.501263), 
             'N24': (0.47163, 0.770632, -0.650495, -0.500592), 
             'N25': (-0.615548, -0.316352, -0.814774, -0.665052), 
             'N26': (-0.305399, -0.00591217, -0.814862, -0.665489), 
             'N27': (0.00550714, 0.304979, -0.815022, -0.665418), 
             'N28': (0.316126, 0.615276, -0.814707, -0.664908), 
             'N29': (-0.46018, -0.16101, -0.97887, -0.829315), 
             'N31': (0.160884, 0.460147, -0.978775, -0.829426), 
             'S1': (-1.08096, -0.782554, 0.00715956, 0.15689), 
             'S2': (-0.7713, -0.47242, 0.0074194, 0.157269), 
             'S3': (-0.4611, -0.161377, 0.00723009, 0.157192), 
             'S4': (-0.149836, 0.150222, 0.00737069, 0.157441), 
             'S5': (0.161297, 0.461031, 0.0072399, 0.1572), 
             'S6': (0.472537, 0.771441, 0.00728934, 0.157137), 
             'S7': (0.782516, 1.08097, 0.00742809, 0.15709), 
             'S8': (-0.92583, -0.627259, 0.171786, 0.32173), 
             'S9': (-0.616329, -0.31694, 0.171889, 0.321823), 
             'S10': (-0.305695, -0.00579187, 0.172216, 0.322179), 
             'S11': (0.00556739, 0.305472, 0.172237, 0.322278), 
             'S12': (0.316973, 0.61631, 0.172015, 0.322057), 
             'S13': (0.627389, 0.925972, 0.171749, 0.321672), 
             'S14': (-0.925847, -0.627123, 0.335898, 0.48578), 
             'S15': (-0.616201, -0.316839, 0.336498, 0.486438), 
             'S16': (-0.305558, -0.00574858, 0.336904, 0.486749), 
             'S17': (0.00557115, 0.305423, 0.33675, 0.486491), 
             'S18': (0.316635, 0.615931, 0.33649, 0.486573), 
             'S19': (0.627207, 0.925969, 0.336118, 0.485923), 
             'S20': (-0.770675, -0.471718, 0.500411, 0.65042), 
             'S21': (-0.46072, -0.161101, 0.501198, 0.650786), 
             'S22': (-0.149915, 0.14982, 0.501334, 0.650856), 
             'S23': (0.160973, 0.460482, 0.501075, 0.650896), 
             'S24': (0.47167, 0.770647, 0.50045, 0.650441), 
             'S25': (-0.615564, -0.316325, 0.66501, 0.814674), 
             'S26': (-0.30512, -0.0056517, 0.665531, 0.81505), 
             'S27': (0.00560886, 0.305082, 0.665509, 0.815022), 
             'S28': (0.316158, 0.615391, 0.665058, 0.814732), 
             'S29': (-0.46021, -0.160988, 0.829248, 0.978699), 
             'S30': (-0.150043, 0.149464, 0.829007, 0.978648), 
             'S31': (0.160898, 0.460111, 0.82932, 0.978804) }

""" Now build another dictionary that gives the CCD centers on the sky in degree system again.
These will be (x,y) tuples.
"""
ccdCenters = {}
for detpos,bounds in ccdBounds.items():
    ccdCenters[detpos] = ( 0.5*(bounds[0]+bounds[1]), 0.5*(bounds[2]+bounds[3]))

def minimalHeader(detpos):
    """ Return a minimal FITS Header object containing the detpos, ccdnum, and the WCS info
    needed to show the mosaic on DS9, plus DETSEC info that give IRAF-style mosaic display.
    """
    ccdnum = ccdnums[detpos]
    h = pf.PrimaryHDU().header
    h['CCDNUM'] = ccdnum
    h['DETPOS'] = detpos
    h['EXTNAME'] = detpos
    h['CRVAL1'] = 0.
    h['CRVAL2'] = 0.
    h['CD1_1'] = 0.
    h['CD1_2'] = +7.286e-5
    h['CD2_2'] = 0.
    h['CD2_1'] = -7.286e-5
    h['CTYPE1'] = 'RA---TAN'
    h['CTYPE2'] = 'DEC--TAN'
    x,y = ccdCorners[detpos]
    h['DETSEC']='[{:d}:{:d},{:d}:{:d}]'.format(x, x+2047, y, y+4095)
    h['CRPIX2']=14826.- ((y-1)//2048)*2129.6667
    h['CRPIX1']=13423.2- ((x-1)//2048)*2254.4
    
    return h
    
class XtalkRemover:
    """Class that can operate on an image HDU from DECAM to remove the intra-CCD
    crosstalk.  Inter-CCD crosstalk is at least 10x smaller and will be ignored for now.
    Default construction will use a built-in set of coefficients.  
    Constructor can also take a string as argument which is name of
    a file in DESDM's xtalk format, from which the coefficients will be moved.
    The corrections will be applied sequentially, and we'll ignore the
    <1e-6 error that results so this can be done in-place.
    """
    def __init__(self, filename=None):
        self.atob = np.zeros(63, dtype=np.float32)
        # Array giving coefficient for A amp source and B amp victim
        self.btoa = np.zeros(63, dtype=np.float32)
        # Array giving coefficient for B amp source and A amp victim

        if filename==None:
            # Here are values from DECam_20130606.xtalk file
            self.atob[ 1]=-6.470e-04; self.btoa[ 1]=+4.240e-04;
            self.atob[ 2]=+1.210e-04; self.btoa[ 2]=+3.120e-04;
            self.atob[ 3]=-6.140e-04; self.btoa[ 3]=+1.570e-04;
            self.atob[ 4]=-4.760e-04; self.btoa[ 4]=+1.480e-04;
            self.atob[ 5]=+2.930e-04; self.btoa[ 5]=+3.400e-04;
            self.atob[ 6]=-4.570e-04; self.btoa[ 6]=+3.820e-04;
            self.atob[ 7]=-9.670e-04; self.btoa[ 7]=+2.820e-04;
            self.atob[ 8]=+1.100e-04; self.btoa[ 8]=+2.500e-04;
            self.atob[ 9]=-8.030e-04; self.btoa[ 9]=+3.560e-04;
            self.atob[10]=-6.130e-04; self.btoa[10]=+3.810e-04;
            self.atob[11]=+9.750e-05; self.btoa[11]=+2.340e-04;
            self.atob[12]=-5.160e-04; self.btoa[12]=+1.690e-04;
            self.atob[13]=-7.520e-04; self.btoa[13]=+3.680e-04;
            self.atob[14]=-6.690e-04; self.btoa[14]=+5.420e-04;
            self.atob[15]=+2.260e-04; self.btoa[15]=+3.520e-04;
            self.atob[16]=-4.930e-04; self.btoa[16]=+4.120e-04;
            self.atob[17]=+1.910e-04; self.btoa[17]=+4.080e-04;
            self.atob[18]=-6.460e-04; self.btoa[18]=+1.270e-04;
            self.atob[19]=+1.560e-04; self.btoa[19]=+3.070e-04;
            self.atob[20]=-5.550e-04; self.btoa[20]=+4.410e-04;
            self.atob[21]=+2.780e-04; self.btoa[21]=+3.860e-04;
            self.atob[22]=-4.920e-04; self.btoa[22]=+3.640e-04;
            self.atob[23]=-9.180e-04; self.btoa[23]=+3.290e-04;
            self.atob[24]=-6.800e-04; self.btoa[24]=+5.140e-04;
            self.atob[25]=+1.340e-04; self.btoa[25]=+3.160e-04;
            self.atob[26]=-5.360e-04; self.btoa[26]=+4.500e-04;
            self.atob[27]=-7.180e-04; self.btoa[27]=+4.470e-04;
            self.atob[28]=+2.330e-04; self.btoa[28]=+4.090e-04;
            self.atob[29]=-3.780e-04; self.btoa[29]=+5.110e-04;
            self.atob[30]=-7.200e-04; self.btoa[30]=+3.160e-04;
            self.atob[31]=-6.710e-04; self.btoa[31]=+4.410e-04;
            self.atob[32]=-8.280e-04; self.btoa[32]=+3.400e-04;
            self.atob[33]=+1.140e-04; self.btoa[33]=+3.280e-04;
            self.atob[34]=-5.230e-04; self.btoa[34]=+1.390e-04;
            self.atob[35]=+2.750e-04; self.btoa[35]=+3.160e-04;
            self.atob[36]=-6.900e-04; self.btoa[36]=+5.530e-04;
            self.atob[37]=-9.650e-04; self.btoa[37]=+4.510e-04;
            self.atob[38]=-8.390e-04; self.btoa[38]=+4.390e-04;
            self.atob[39]=+3.410e-04; self.btoa[39]=+6.580e-04;
            self.atob[40]=-5.700e-04; self.btoa[40]=+4.630e-04;
            self.atob[41]=+1.520e-04; self.btoa[41]=+4.440e-04;
            self.atob[42]=-6.040e-04; self.btoa[42]=+4.830e-04;
            self.atob[43]=-8.870e-04; self.btoa[43]=+3.760e-04;
            self.atob[44]=-6.560e-04; self.btoa[44]=+4.050e-04;
            self.atob[45]=-7.120e-04; self.btoa[45]=+4.840e-04;
            self.atob[46]=-6.840e-04; self.btoa[46]=+4.540e-04;
            self.atob[47]=+9.290e-05; self.btoa[47]=+4.270e-04;
            self.atob[48]=-6.090e-04; self.btoa[48]=+5.570e-04;
            self.atob[49]=+4.900e-05; self.btoa[49]=+3.570e-04;
            self.atob[50]=-6.280e-04; self.btoa[50]=+1.070e-04;
            self.atob[51]=+1.120e-04; self.btoa[51]=+3.400e-04;
            self.atob[52]=-9.240e-04; self.btoa[52]=+4.780e-04;
            self.atob[53]=-7.890e-04; self.btoa[53]=+4.990e-04;
            self.atob[54]=+5.960e-05; self.btoa[54]=+3.820e-04;
            self.atob[55]=-7.430e-04; self.btoa[55]=+0.000e+00;  # missing value here?
            self.atob[56]=-4.880e-04; self.btoa[56]=+1.730e-04;
            self.atob[57]=+1.160e-04; self.btoa[57]=+4.370e-04;
            self.atob[58]=-6.660e-04; self.btoa[58]=+5.060e-04;
            self.atob[59]=-9.930e-04; self.btoa[59]=+3.350e-04;
            self.atob[60]=-6.860e-04; self.btoa[60]=+4.480e-04;
            self.atob[61]=+1.120e-04; self.btoa[61]=+3.330e-04;
            self.atob[62]=-7.120e-04; self.btoa[62]=+1.040e-04;

        else:
            # Read the values from file
            f = open(filename)
            for line in f:
                if len(line.strip())==0 or line[0]=='#': continue
                fields = line.split()
                m = re.match(r"ccd(\d*)([AB])",fields[0])
                if not m:
                    raise Exception('XtalkRemover did not get good victim field: ' +fields[0])
                ccdnum1 = int(m.group(1))
                amp1 = m.group(2)

                m = re.match(r"ccd(\d*)([AB])",fields[1])
                if not m:
                    raise Exception('XtalkRemover did not get good source field: ' + fields[0])
                ccdnum2 = int(m.group(1))
                amp2 = m.group(2)

                if ccdnum1==ccdnum2:
                    # Use intra-chip xtalk coefficients
                    if amp1=='A' and amp2=='B':
                        self.btoa[ccdnum1] = float(fields[2])
                    elif amp1=='B' and amp2=='A':
                        self.atob[ccdnum1] = float(fields[2])
                    else:
                        raise Exception('XtalkRemover got bad amp specs: '+amp1+' -> '+amp2)
        return

    def __call__(self, hdu):
        """Use as an operator on HDU of an image from DECAM to apply the xtalk correction
        """
        detpos = hdu.header['DETPOS']
        ccdnum = ccdnums[detpos]
        shape = hdu.data.shape
        if len(shape) != 2:
            raise Exception('XtalkRemover() got hdu that is not 2d image')
        xCenter = shape[1]//2
        if detpos[0]=='N':
            # North CCDs have A at lower x
            ltor = self.atob[ccdnum]
            rtol = self.btoa[ccdnum]
        elif detpos[0]=='S':
            # South CCDs have B at lower x values
            ltor = self.btoa[ccdnum]
            rtol = self.atob[ccdnum]
                
        hdu.data[:,:xCenter] -= hdu.data[:, :xCenter-1:-1] * rtol
        hdu.data[:,xCenter:] -= hdu.data[:, xCenter-1::-1] * ltor
        return
