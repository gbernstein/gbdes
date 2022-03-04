from __future__ import division,print_function
import numpy as np
import numpy.ma as ma

def clippedMean(a, nSigma, axis=None, sigma=None, sigmaFloor=None, maxSample=None):
    """ Return (mean, variance, n) of a vector after percentile-based sigma clip
    sigma is defined from inter-quartile range.  Points outside +-nSigma x IQD/1.349
    are removed before calculating mean, variance, and number of points.
    If sigma is given, clip outside nSigma x sigma from the median
    If sigmaFloor is given, it sets a lower limit to clipping sigma
    a      = array of data
    nSigma = number of sigma at which to clip
    axis   = axis along which all stats to be taken.  If None, uses entire flattened array
    sigma  = sigma of the data.  If None, determines from inter-quartile range
    sigmaFloor = lowest allowable sigma
    maxSample  = maximum number of points to use in finding clip limits.  If None, use all.

    Returns:
    mean   = mean of unclipped points
    var    = variance of unclipped points
    n      = number of unclipped points.
    Each is scalar if axis==None, or 1 lower dimension than input array.
    """
    if axis==None:
        dataLength = len(a.flatten())
    else:
        dataLength = a.shape[axis]
    if maxSample==None or dataLength <= maxSample:
        sample = a
    else:
        # Subsample the data in deriving clipping limits
        step = (dataLength-1)//maxSample + 1
        if axis==None:
            sample = a.flatten[::step]
        else:
            s = []
            for i in range(0,axis):
                s.append(slice(None))
            s.append(slice(None,None,step))
            for i in range(axis+1,len(a.shape)):
                s.append(slice(None))
            sample = a[s]
                
    if sigma==None:
        # Determine sigma from IQD
        iqdSigma = 1.349
        p75 = np.percentile(sample, 75., axis=axis)
        p25 = np.percentile(sample, 25., axis=axis)
        if sigmaFloor==None:
            # Use the IQD straight up:
            upper = (0.5+nSigma/iqdSigma)*p75 + (0.5-nSigma/iqdSigma)*p25
            lower = (0.5-nSigma/iqdSigma)*p75 + (0.5+nSigma/iqdSigma)*p25
        else:
            # set sigma as maximum of IQD and the floor
            dev = nSigma*np.maximum( (p75-p25)/iqdSigma, sigmaFloor)
            mid = 0.5*(p25+p75)
            lower = mid-dev
            upper = mid+dev
    else:
        # use prescribed sigma about median
        med = np.median(sample, axis=axis)
        upper = med + nSigma*sigma
        lower = med - nSigma*sigma

    if axis==None:
        data = ma.masked_outside(a, lower, upper)
    else:
        # Need to broadcast the upper and lower bounds across chosen axis
        bshape = a.shape[:axis] + (1,) + a.shape[axis+1:]
        mask = np.logical_or(a < lower.reshape(bshape), a>upper.reshape(bshape))
        data = ma.masked_where(mask, a, copy=False)
    if axis==None:
        # Force return of double precision to use this for accumulating mean & avoid roundoff
        return data.mean(dtype=np.float64), data.var(dtype=np.float64), data.count()
    else:
        # A problem is that count() appears to return float data if it returns an array.
        # Also need to force calculation of the mean and variance in double precision
        return data.mean(axis=axis, dtype=np.float64).astype(np.float32), \
          data.var(axis=axis, dtype=np.float64).astype(np.float32), \
          data.count(axis=axis).astype(np.int32)
