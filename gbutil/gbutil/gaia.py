# From Alex Drlica Wagner
# https://cdcvs.fnal.gov/redmine/projects/des-sci-release/repository/entry/users/kadrlica/catalog_coadd/trunk/code/utils.py#L275
import warnings
#from astroquery.vizier import Vizier
from astropy.coordinates import SkyCoord
import astropy.units as u
import astropy.io.fits as pf
import numpy as np

def get_vizier_catalog(ra,dec,radius=None,**kwargs):

    kwargs.setdefault('row_limit',-1)
    coord = SkyCoord(ra*u.deg,dec*u.deg)
    radius = u.Quantity(radius,u.deg)
    vizier = Vizier(**kwargs)
    warnings.filterwarnings("ignore")
    tab = vizier.query_region(coord,radius)
    warnings.resetwarnings()
    return tab[0]


def gaia_to_fits(ra, dec, radius, fitsname='gaia.cat'):
    gaia = dict(catalog='I/337/gaia',
                columns=['RA_ICRS','DE_ICRS','e_RA_ICRS','e_DE_ICRS','<Gmag>'])
    c = get_vizier_catalog(ra,dec,radius,**gaia)
    err = np.maximum(c['e_RA_ICRS'],c['e_DE_ICRS'])
    maximum_gaia_error = 5.  # Largest uncertainty to keep
    use = np.where(err < 5.)
    # Convert errors from mas to degrees for table
    hdu = pf.BinTableHDU.from_columns( [ \
        pf.Column(name='RA_ICRS', format='D', array=c['RA_ICRS'][use]),
        pf.Column(name='DE_ICRS', format='D', array=c['DE_ICRS'][use]),
        pf.Column(name='ERROR', format='E', array=err[use] * 0.001 / 3600.),
        pf.Column(name='GMAG', format='E',array=c['__Gmag_'][use])])
    hdu.header['RA'] = ra
    hdu.header['Dec'] = dec
    hdu.header['RADIUS'] = radius
    hdu.writeto(fitsname, overwrite=True)
    return
