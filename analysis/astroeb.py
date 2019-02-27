# Routines for E/B separations of astrometric residuals
# and investigation of correlation functions etc.
# Note the use of Mike Jarvis's treecorr below
from __future__ import division,print_function
import numpy as np
import astropy.io.fits as pf
import astropy.wcs as wcs
import pylab as pl
import re
import treecorr as tc
import astroplots as ap

def vcorr(x,y,dx,dy,
          rmin=5./3600., rmax=1.5, dlogr=0.05,
          maxpts = 30000):
    """
    Produce angle-averaged 2-point correlation functions of astrometric error
    for the supplied sample of data, using brute-force pair counting.
    Output are the following functions:
    logr - mean log of radius in each bin
    xi_+ - <vr1 vr2 + vt1 vt2> = <vx1 vx2 + vy1 vy2>
    xi_- - <vr1 vr2 - vt1 vt2>
    xi_x - <vr1 vt2 + vt1 vr2>
    xi_z2 - <vx1 vx2 - vy1 vy2 + 2 i vx1 vy2>
    """

    if len(x) > maxpts:
        # Subsample array to get desired number of points
        rate = float(maxpts) / len(x)
        print("Subsampling rate {:5.3f}%".format(rate*100.))
        use = np.random.random(len(x)) <= rate
        x = x[use]
        y = y[use]
        dx = dx[use]
        dy = dy[use]
    print("Length ",len(x))
    # Get index arrays that make all unique pairs
    i1, i2 = np.triu_indices(len(x))
    # Omit self-pairs
    use = i1!=i2
    i1 = i1[use]
    i2 = i2[use]
    del use
    
    # Make complex separation vector
    dr = 1j * (y[i2]-y[i1])
    dr += x[i2]-x[i1]

    # log radius vector used to bin data
    logdr = np.log(np.absolute(dr))
    logrmin = np.log(rmin)
    bins = int(np.ceil(np.log(rmax/rmin)/dlogr))
    hrange = (logrmin, logrmin+bins*dlogr)
    counts = np.histogram(logdr, bins=bins, range=hrange)[0]
    logr = np.histogram(logdr, bins=bins, range=hrange, weights=logdr)[0] / counts

    # First accumulate un-rotated stats
    v =  dx + 1j*dy
    vv = dx[i1] * dx[i2] + dy[i1] * dy[i2]
    xiplus = np.histogram(logdr, bins=bins, range=hrange, weights=vv)[0]/counts
    vv = v[i1] * v[i2]
    xiz2 = np.histogram(logdr, bins=bins, range=hrange, weights=vv)[0]/counts

    # Now rotate into radial / perp components
    vv *= np.conj(dr)
    vv *= np.conj(dr)
    dr = dr.real*dr.real + dr.imag*dr.imag
    vv /= dr
    del dr
    ximinus = np.histogram(logdr, bins=bins, range=hrange, weights=vv)[0]/counts
    xicross = np.imag(ximinus)
    ximinus = np.real(ximinus)

    return logr, xiplus, ximinus, xicross, xiz2

def xiB(logr, xiplus, ximinus):
    """
    Return estimate of pure B-mode correlation function
    """
    # Integral of d(log r) ximinus(r) from r to infty:
    dlogr = np.zeros_like(logr)
    dlogr[1:-1] = 0.5*(logr[2:] - logr[:-2])
    tmp = np.array(ximinus) * dlogr
    integral = np.cumsum(tmp[::-1])[::-1]
    return 0.5*(xiplus-ximinus) + integral

def ebSums(fitsname, exposures):
    """
    Do brute-force correlations and eb split for a list of exposures.
    """
    ff = pf.open(fitsname)
    outfile = open('ebsums.dat','w')
    sig = 0.02 # sigma of the Gaussian we'll use for 2nd moments of power specs
    nsum = 0
    for exp in exposures:
        t = ap.select(ff, exposure=exp, clip=False)
        x = t['xW']
        y = t['yW']
        dx = t['xresW'] / (1-t['wtfrac'])
        dy = t['yresW'] / (1-t['wtfrac'])
        # Filter low-order poly out of data
        polyfit(x,y,dx,dy,order=3)
        out = vcorr(x,y,dx,dy,maxpts=20000)
        if nsum==0:
            logr = out[0]
            esum = np.zeros_like(logr)
            bsum = np.zeros_like(logr)
            xsum = np.zeros_like(logr)
            r = np.exp(logr)
            wr = r * np.exp(-0.5*r*r/(sig*sig))
            wr /= sum(wr)
        xib = xiB(out[0], out[1], out[2])
        xie = out[1] - xib
        pl.plot(out[0], xie, label=exp)
        pl.plot(out[0], xib)
        esum += xie
        bsum += xib
        xsum += out[3]
        nsum += 1
        var = np.sum(out[1]*wr)
        ee = np.sum(out[4]*wr) / var
        print(exp, var, np.absolute(ee), np.angle(ee,deg=True)/2.)
        print(exp, var, np.absolute(ee), np.angle(ee,deg=True)/2.,file=outfile)
    outfile.close()
    return logr, esum/nsum, bsum/nsum, xsum/nsum

def fastXi(x,y,dx,dy,rrange=(5./3600.,1.5), nbins=100):
    """
    Use treecorr to get xi of (dx*dx + dy*dy).
    Returns arrays of mean radius in each (log-spaced) bin, mean xi.
    """
    catx = tc.Catalog(x=x, y=y, k=dx)
    caty = tc.Catalog(x=x, y=y, k=dy)
    kk = tc.KKCorrelation(min_sep=rrange[0], max_sep=rrange[1], nbins=nbins)
    kk.process(catx)
    xx = np.array(kk.xi)

    kk = tc.KKCorrelation(min_sep=rrange[0], max_sep=rrange[1], nbins=nbins)
    kk.process(caty)
    yy = np.array(kk.xi)

    return np.exp(kk.meanlogr), xx+yy

def vcorr2d(x,y,dx,dy,
            rmax=1., bins=513,
            maxpts = 30000):
    """
    Produce 2d 2-point correlation function of total displacement power
    for the supplied sample of data, using brute-force pair counting.
    Output are 2d arrays giving the 2PCF and then the number of pairs that
    went into each bin.  The 2PCF calculated is
    xi_+ - <vr1 vr2 + vt1 vt2> = <vx1 vx2 + vy1 vy2>
    Note that each pair is counted only once.  So to count all pairs one can
    average xi_+ with itself reflected about the origin.
    """

    hrange = [ [-rmax,rmax], [-rmax,rmax] ]
    if len(x) > maxpts:
        # Subsample array to get desired number of points
        rate = float(maxpts) / len(x)
        print("Subsampling rate {:5.3f}%".format(rate*100.))
        use = np.random.random(len(x)) <= rate
        x = x[use]
        y = y[use]
        dx = dx[use]
        dy = dy[use]
    print("Length ",len(x))
    # Get index arrays that make all unique pairs
    i1, i2 = np.triu_indices(len(x))
    # Omit self-pairs
    use = i1!=i2
    i1 = i1[use]
    i2 = i2[use]
    del use
    
    # Make separation vectors and count pairs
    yshift = y[i2]-y[i1]
    xshift = x[i2]-x[i1]
    counts = np.histogram2d(xshift,yshift, bins=bins, range=hrange)[0]

    # Accumulate displacement sums
    v =  dx + 1j*dy
    print('xiplus') ##
    vv = dx[i1] * dx[i2] + dy[i1] * dy[i2]
    xiplus = np.histogram2d(xshift,yshift, bins=bins, range=hrange, weights=vv)[0]/counts

    return xiplus,counts

def polyfit(x,y,dx,dy,order=3):
    """
    Fit and remove a polynomial function of (x,y) from (dx,dy).
    """
    xpow = np.ones_like(x)
    terms = []
    for i in range(0,order+1):
        term = np.array(xpow);
        terms.append(np.array(term))
        for j in range(1,order+1-i):
            term *= y
            terms.append(np.array(term))
        xpow *= x
    xy = np.transpose(np.vstack(terms))
    dxy = np.transpose(np.vstack( (dx,dy)))
    # Fit and subtract dx, dy concurrently
    a = np.linalg.lstsq(xy, dxy)[0]
    print("Coefficients: ",a)
    fit = np.dot(xy,a)
    dx -= fit[:,0]
    dy -= fit[:,1]
    return

def ebSum2d(fitsname, exposures):
    """
    Do brute-force correlations and sum pairs over many exposures to make a good 2d 2PCF
    """
    ff = pf.open(fitsname)
    nsum = 0
    for exp in exposures:
        t = ap.select(ff, exposure=exp, clip=False)
        x = t['xW']
        y = t['yW']
        dx = t['xresW'] / (1-t['wtfrac'])
        dy = t['yresW'] / (1-t['wtfrac'])
        # Filter low-order poly out of data
        polyfit(x,y,dx,dy,order=3)
        x,c = vcorr2d(x,y,dx,dy,maxpts=20000)
        if nsum==0:
            cxi = np.zeros_like(x)
            counts = np.zeros_like(c)
        cxi += x*c
        counts += c
        nsum += 1
    return cxi/counts, counts

def summary(x,y,dx,dy, order=3, sigma=1./60.):
    """
    Return some summary statistics of the residual vector field given by x,y,dx,dy.
    First fit and subtract the best-fitting polynomial of the selected order.
    Then calculate the xi of the remaining residuals.  Return the correlation function
    value averaged near the origin with Gaussian weight of the selected sigma (1 arcmin by default).
    Then find the radius at which the correlations drop to half of this value.
    Returns:
      xi0:  Mean <dx*dx+dy*dy> for pairs at small separation 
      rCorr: Radius at which correlations drop to half of this
      polyVar: Variance (over the sample points) of the best-fit polynomial
      xcoeffs: Coefficients of the polynomial that was fit and subtracted
      ycoeffs: Coefficients of the polynomial that was fit and subtracted
    """
    xp = ap.Poly2d(order, sumOrder=True)
    yp = ap.Poly2d(order, sumOrder=True)
    xp.fit(x,y,dx)
    xfit = xp.evaluate(x,y)
    yp.fit(x,y,dy)
    yfit = yp.evaluate(x,y)
    polyVar = np.var(xfit) + np.var(yfit)

    r, xi = fastXi(x,y,dx-xfit,dy-yfit)
    wt = r * r * np.exp(-0.5*r*r/sigma/sigma)
    xi0 = np.sum(wt*xi) / np.sum(wt)

    # Find last radius with xi>0.5*xi0:
    lastR = np.where(xi>0.5*xi0)[0][-1]
    # Interpolate the crossing point
    f = (0.5*xi0 - xi[lastR]) / (xi[lastR+1]-xi[lastR])
    rCorr = r[lastR] * (r[lastR+1]/r[lastR])**f
    return xi0, rCorr, polyVar, xp.getCoeffs(), yp.getCoeffs()

def allExposuresSummary(fitsname, savefile, order=3):
    """
    Collect astrometric error summary information from all exposures in the wcsfit
    output file with name fitsname.  Save in a binary fits table in file savefile.
    """
    ff = pf.open(fitsname)
    # Get all of our exposures
    exptab = ff['EXPOSURES'].data

    # Make a dictionary of band names per instrument
    bands = {}
    for hdu in ff[1:]:
        if hdu.header['EXTNAME']=='Instrument':
            bands[hdu.header['NUMBER']] = hdu.header['BAND']

    nexp = np.count_nonzero(exptab['INSTRUMENTNUMBER']>=0)

    # Set up the table we want to build
    coeff_format = '{:d}E'.format((order+1)*(order+2)//2) # FITS format code for the array of poly coeffs
    hdu = pf.BinTableHDU.from_columns( \
        (pf.Column(name='exposure',format='20A'),
         pf.Column(name='band',format='2A'),
         pf.Column(name='MJD',format='D'),
         pf.Column(name='exptime',format='D'),
         pf.Column(name='xi0',format='E'),
         pf.Column(name='rCorr',format='E'),
         pf.Column(name='polyvar',format='E'),
         pf.Column(name='xcoeff',format=coeff_format),
         pf.Column(name='ycoeff',format=coeff_format)),
        nrows = nexp)
    
    # ap.select is too slow to do every exposure so I'm going to speed it up knowing
    # we want to use basically all objects anyway
    objtab = ff['WCSOut'].data
    print("..read",len(objtab)) ##
    use = objtab['CLIP']==0
    use = np.logical_and(use, objtab['WTFRAC'] < 0.35)
    objtab = objtab[use]
    print("..using",len(objtab)) ##
    # Sort the table by increasing extension number
    ii = np.argsort(objtab['EXTENSION'])
    objtab = objtab[ii]
    print("..sorted") ##
    # Locations of first ii for each extension
    inew = np.where(objtab['EXTENSION'][1:]!=objtab['EXTENSION'][:-1])[0] + 1  
    print("..len(inew)",len(inew)) ##
    extnSlices = {}
    for j in range(len(inew)-1):
        extnSlices[objtab['EXTENSION'][inew[j]]] = slice(inew[j],inew[j+1])
    extnSlices[objtab['EXTENSION'][inew[-1]]] = slice(inew[-1])  # Last extension
    del inew
    print(len(extnSlices), 'slices')

    # Put each exposure's data into the table
    # Process the exposures in order of their names
    names = np.array(["".join(n) for n in exptab['NAME']])
    #ii = np.argsort(exptab['NAME'])
    ii = np.argsort(names)
    iout = -1
    for i in ii:
        if exptab['INSTRUMENTNUMBER'][i]<0:
            continue
        iout = iout + 1
        #expname = exptab['NAME'][i]
        expname = names[i]
        mjd = exptab['MJD'][i]
        exptime = exptab['EXPTIME'][i]
        band = bands[exptab['INSTRUMENTNUMBER'][i]]

        ##s = ap.select(ff, exposure=expname)
        ##out = summary(s['xw'],s['yw'],s['xresw'],s['yresw'],order=order)
        # Get data from all extensions using this exposure
        extns = np.where(ff['EXTENSIONS'].data['EXPOSURE']==i)[0]
        xx = []
        yy = []
        dxx = []
        dyy = []
        for extn in extns:
            if extn not in extnSlices:
                continue
            xx.append(objtab['xw'][extnSlices[extn]])
            yy.append(objtab['yw'][extnSlices[extn]])
            dxx.append(objtab['xresw'][extnSlices[extn]])
            dyy.append(objtab['yresw'][extnSlices[extn]])
        if len(xx) < 10:
            print("Skipping",expname,"with too few matched extensions: ",len(xx))
            continue
        x = np.concatenate(xx)
        if len(x) < 2000:
            print("Skipping",expname,"with too few stars: ",len(x))
        y = np.concatenate(yy)
        dx = np.concatenate(dxx)
        dy = np.concatenate(dyy)
        
        out = summary(x,y,dx,dy, order=order)
        hdu.data[iout] = (expname, band, mjd, exptime) + out  # Write whole table line
        print(expname, out[:3])

    hdu.writeto(savefile, clobber=True)
    return

    
    
