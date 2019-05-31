import astropy.io.fits as pf
import astropy.table as tb
import numpy as np 
import treecorr as tc 
import astroplots as ap

def get_Expo(fits):
    '''
    This function finds all the exposures in a given SAVE residual file, coming from Gary's astrometric codes. This is the first step into computing xi
    '''

    if type(fits)==str:
        ff = pf.open(fits)
    else:
        ff = fits
    #As the exposures are saved in a [D, X, X, X, X , X, X, X] format, we should have a way of copying this...    
    a = [''.join(i['Name']) for i in ff['Exposures'].data] #magical one line double for 
    return a

def get_data(fits, expo):
    '''
    Copies the data needed for both xi and sigma
    '''
    res = ap.select(fits,exposure=expo)
    cov = res['covTotalW']
    wt  = 2. / np.sum(cov[:,:2]*cov[:,:2],axis=1)
    wtAdjust = 2. / res['chisqExpected']
    xv = res['xresW'] * wtAdjust
    yv = res['yresW'] * wtAdjust
    xW = res['xW']
    yW = res['yW']

    return xv, yv, wt, xW, yW

def data_dump(fits, exposure):
    ''' This copies the residuals from the SAVE file and saves in a manner that xi_catalog_tc can use '''
    xv, yv, sig, xW, yW = get_data(fits, exposure)
    a = tb.Table()
    a['DELTA_X'] = xv
    a['DELTA_Y'] = yv
    a['X'] = xW
    a['Y'] = yW

    return a

def xi_catalog_tc(fits, output, num_threads, return_data = False):
    '''
    Builds a catalog of <xi_x>, <xi_y> and <xi_c> for a given set of exposures using treecorr
    '''
    fits = pf.open(fits)
    exposures = get_Expo(fits)
    table = tb.Table(names=("EXPNUM", "XI_X", "XI_Y", "XI_C", 'N'))
    for i in exposures:
        try:
            data = data_dump(fits, i, return_data = True)
            av_x, av_y, av_c = compute_xi(data['X'], data['Y'], data['DELTA_X'], data['DELTA_Y'], num_threads)
            table.add_row([int(i[1:]), av_x, av_y, av_c, len(data)])
            del data
        except:
            exposures.remove(i)

    table.sort('EXPNUM')
    table.write(output)

    if return_data:
        return table


def compute_xi(x, y, delta_x, delta_y, num_threads):
    x_cat = tc.Catalog(x = x, y = y, k = delta_x)
    y_cat = tc.Catalog(x = x, y = y, k = delta_y)
    kk = tc.KKCorrelation(min_sep=5./3600, max_sep=1.5, nbins=200)
    kk.process(x_cat, num_threads = num_threads)
    xi_x, logr = kk.xi, kk.logr         
    n = kk.npairs[(logr > -5) & (logr < -4.5)]
    av_x =  np.average(xi_x[(logr > -5) & (logr < -4.5)], weights=n)
    kk.process(y_cat, num_threads = num_threads)
    xi_y, logr = kk.xi, kk.logr         
    n = kk.npairs[(logr > -5) & (logr < -4.5)]
    av_y =  np.average(xi_y[(logr > -5) & (logr < -4.5)], weights=n)
    kk.process(x_cat, y_cat, num_threads = num_threads)
    xi_c, logr = kk.xi, kk.logr         
    n = kk.npairs[(logr > -5) & (logr < -4.5)]
    av_c =  np.average(xi_c[(logr > -5) & (logr < -4.5)], weights=n)

    return av_x, av_y, av_c


if __name__ == '__main__':
    import sys

    filename = sys.argv[1]
    output = sys.argv[2]
    try:
        num_threads = int(sys.argv[3])
    except:
        num_threads = 8


    xi_catalog_tc(filename, output, num_threads)


