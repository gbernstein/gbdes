'''
Implementation of gbtools/astrometry PixelMap classes in Python.
This is meant as a read-only version of the astrometry library.
Capabilities include:
* deserialize the YAML files specifying PixelMaps and WCS's
* Execute forward (pixel to world) transformations with any map/wcs
* Calculate local derivatives of maps
* Execute inverse transformations using generic solver method

All transformations can be done on arrays of coordinates as well as scalars.

Pixel or world positions are to be numpy arrays of shape (2) or (N,2).  Sky
positions are astropy.coordinates.SkyCoord objects/arrays.
'''

import numpy as np
import astropy.coordinates as co
import future
import yaml
from scipy.optimize import root


class PixelMap(object):
    ''' Base class for transformations from one 2d system ("pixel") to another ("world").
    Each derived class must implement __call__ to execute this transform on array of
    shape (2) or (N,2).  Base class implements some routines such as taking local derivatives
    and solving for map inverse.
    c is always a color (or array of them)
    '''
    def __init__(self, name):
        self.name = name
        return
    def jacobian(self, xy, c=None, step=1.):
        '''Use finite differences with indicated step size to get Jacobian derivative
        matrix at each point in xy.  If xy.shape=(N,2) then the returned array has
        shape (N,2,2) with [n,i,j] giving dx_i/dx_j at nth point.
        '''
        dx = (self(xy+np.array([step,0.])) - self(xy-np.array([step,0.]))) / (2*step)
        dy = (self(xy+np.array([0.,step])) - self(xy-np.array([0.,step]))) / (2*step)
        npts = np.atleast_2d(xy).shape[0]
        out = np.zeros((npts,2,2),dtype=float)
        out[:,:,0] = dx
        out[:,:,1] = dy
        return out
    def invert(self, xyw, xyp, c=None, tol=0.0001):
        '''
        Fill the array xyp with the solutions to PixelMap(xyp)=xyw.  The input values
        of xyp are the initial guess.  Tolerance can be altered from default, which
        is given in the pixel space ???
        '''
        # Need to call the solver row by row, will be slow.
        class resid(object):
            ''' Callable giving deviation of output from target
            '''
            def __init__(self,map,target,c=None):
                self.target = target
                self.map = map
                self.c = c
                return
            def __call__(self,xyp):
                return self.map(xyp,self.c) - self.target
        xyw2d = np.atleast_2d(xyw)
        xyp2d = np.atleast_2d(xyp)
        npts = xyw2d.shape[0]
        if xyp2d.shape[0] != npts:
            raise Exception('Mismatch of xyw, xyp point counts in PixelMap.invert')

        for i in range(npts):
            result = root(resid(self.__call__,xyw2d[i],c), xyp2d[i], tol=tol).x
            xyp2d[i] = result
        return
    
'''
Following are various derived classes of PixelMaps.  They each have a class attribute
field "type" which corresponds to the "Type" field in the YAML node.
'''
    
class Identity(PixelMap):
    '''Identity map
    '''
    @staticmethod
    def type():
        return 'Identity'
    def __init__(self,name, **kwargs):
        super(Identity,self).__init__(name)
        return
    def __call__(self,xy,c=None):
        return np.array(xy)

class Constant(PixelMap):
    '''Constant shift pixel map
    '''
    @staticmethod
    def type():
        return 'Constant'
    def __init__(self, name, **kwargs):
        super(Constant,self).__init__(name)
        if 'Parameters' not in kwargs:
            raise Exception('Missing Parameters in Constant PixelMap spec')
        if len(kwargs['Parameters']) !=2:
            raise Exception('Wrong # of parameters in Constant PixelMap spec')
        self.shift = np.array(kwargs['Parameters'],dtype=float)
        return
    def __call__(self,xy,c=None):
        return xy + self.shift

class Linear(PixelMap):
    '''Affine transformation pixel map
    '''
    @staticmethod
    def type():
        return 'Linear'
    def __init__(self, name, **kwargs):
        super(Constant,self).__init__(name)
        if 'Parameters' not in kwargs:
            raise Exception('Missing Parameters in Linear PixelMap spec')
        if len(kwargs['Parameters']) !=6:
            raise Exception('Wrong # of parameters in Linear PixelMap spec')
        p = np.array(kwargs['Parameters'],dtype=float).reshape((2,3))
        self.shift = p[:,0]
        self.mT = p[:,1:].T  # Make transpose of the magnification matrix
        return
    def __call__(self,xy,c=None):
        return np.dot(xy,self.mT) + self.shift

class Polynomial(PixelMap):
    '''2d polynomial pixel map
    '''
    @staticmethod
    def type():
        return 'Poly'
    def __init__(self,name,**kwargs):
        super(Polynomial,self).__init__(name)
        # Read the scaling bounds that will map interval into [-1,+1]
        xmin = kwargs['XMin']
        xmax = kwargs['XMax']        
        ymin = kwargs['YMin']
        ymax = kwargs['YMax']        
        self.shift = np.array( [0.5*(xmin+xmax), 0.5*(ymin+ymax)])
        self.scale = np.array( [2./(xmax-xmin), 2./(ymax-ymin)])
        # Read coefficients, x then y
        self.coeffs = []
        for d in (kwargs['XPoly'], kwargs['YPoly']):
            sumOrder = bool(kwargs['SumOrder'])
            xOrder = int(kwargs['XOrder'])
            if sumOrder:
                yOrder = xOrder
            else:
                yOrder = int(kwargs['YOrder'])
            c = np.zeros( (xOrder+1, yOrder+1), dtype=float)
            v = np.array(kwargs['Coefficients'])
            # Here is where we are very specific to the ordering of coeffs
            # in the gbtools/utilities/Poly2d.cpp code:
            if sumOrder:
                # Ordering is 1, x, y, x^2, xy, y^2, x^3, xy^2
                ix = 0
                iy = 0
                N = 0
                for cc in v:
                    c[ix,iy] = cc
                    if ix>0:
                        iy = iy+1
                        ix = ix-1
                    else:
                        N = N+1
                        ix = N
                        iy = 0
            else:
                # ordering is 1., y, y^2, y^3 ... , x, xy, xy^2,
                c = v.reshape(xOrder, yOrder)
            self.coeffs.append(np.array(c))
        return
    def __call__(self, xy, c=None):
        xyscaled = (xy - self.shift)*self.scale
        xyw = np.zeros_like(xy)
        xyw[:,0] = np.polynomial.polynomial.polyval2d(xyscaled[:,0],xyscaled[:,1],
                                                     self.coeffs[0])
        xyw[:,1] = np.polynomial.polynomial.polyval2d(xyscaled[:,0],xyscaled[:,1],
                                                     self.coeffs[1])
        return xyw
        
class Template(PixelMap):
    '''Mapping defined by scaling of a 1-dimensional lookup table
    '''
    @staticmethod
    def type():
        return 'Template'

    # A class variable contains the templates read from the last reference file.
    libraryFile = None
    library = {}
    
    def __init__(self,name, **kwargs):
        super(Piecewise,self).__init__(name)
        if kwargs['HasSplit']:
            raise Exception('Template pixel map not coded for split templates')
        fname = kwargs['Filename']
        if libraryFile is None or libraryFile!=fname:
            # Read in a new library file
            libraryFile = fname
            library = yaml.load(open(fname))

        # Now find the desired template
        if kwargs['LowTable'] not in library:
            raise Exception('Did not find map ' + kwargs['LowTable'] + ' in file ' + fname)
        self.scale = float(kwargs['Parameter'])
        tab = library[kwargs['LowTable']]
        
        self.axis = tab['Axis']  # 'X', 'Y', 'R"
        if self.axis=='R':
            self.center = np.array([tab['XCenter'],tab['YCenter']], dtype=float)
        self.vals = np.array(tab['Parameters'])
        # First and last elements are always set to zero in C++:
        self.vals[0] = 0.  ### ??? Check bounds for Template maps !!
        self.vals[-1] = 0.
        self.args = float(tab['ArgStart']) + float(tab['ArgStep']) * \
          np.arange(self.vals.shape[0])
        return

    def __call__(self, xy, c=None):
        # Note that np.interp does not exploit the equal spacing of arguments???
        xyw = np.array(np.atleast_2d(xy))
        if self.axis=='X':
            xyw[:,0] += np.interp(xyw[:,0], self.args, self.vals) * self.scale
        elif self.axis=='Y':
            xyw[:,1] += np.interp(xyw[:,1], self.args, self.vals) * self.scale
        elif self.axis=='R':
            xyc = xyw - self.center
            rad = np.hypot(xyc[:,0],xyc[:,1])
            dr = np.interp(rad, self.args, self.vals) * self.scale
            xyw = xyc*(dr/rad) + xy
        else:
            raise Exception('Unknown Template axis type ' + self.axis)
        return xyw

class Piecewise(PixelMap):
    '''Mapping defined by piecewise-linear function
    '''
    @staticmethod
    def type():
        return 'Piecewise'
    def __init__(self,name,**kwargs):
        super(Piecewise,self).__init__(name)
        self.axis = kwargs['Axis']  # 'X', 'Y', 'R"
        if self.axis=='R':
            self.center = np.array([kwargs['XCenter'],kwargs['YCenter']], dtype=float)
        self.vals = np.array(kwargs['Parameters'])
        # First and last elements are always set to zero in C++:
        self.vals[0] = 0.
        self.vals[-1] = 0.
        self.args = float(kwargs['ArgStart']) + float(kwargs['ArgStep']) * \
          np.arange(self.vals.shape[0])
    def __call__(self, xy, c=None):
        # Note that np.interp does not exploit the equal spacing of arguments???
        xyw = np.array(np.atleast_2d(xy))
        if self.axis=='X':
            xyw[:,0] += np.interp(xyw[:,0], self.args, self.vals)
        elif self.axis=='Y':
            xyw[:,1] += np.interp(xyw[:,1], self.args, self.vals)
        elif self.axis=='R':
            xyc = xyw - self.center
            rad = np.hypot(xyc[:,0],xyc[:,1])
            dr = np.interp(rad, self.args, self.vals)
            xyw = xyc*(dr/rad) + xy
        else:
            raise Exception('Unknown Piecewise axis type ' + self.axis)
        return xyw
        
    
class ColorTerm(PixelMap):
    '''Pixel map that is color times deviation of some other PixelMap
    '''
    @staticmethod
    def type():
        return 'Color'
    def __init__(self, name, map=None, **kwargs):
        '''Need to provide the modified atomic PixelMap as map argument.
        '''
        super(ColorMap,self).__init__(name)
        if map is None:
            raise Exception('No modified map specified for ColorTerm')
        self.map = map
        self.reference = float(kwargs['Reference'])
    def __call__(self, xy, c=None):
        if c is None:
            raise Exception('ColorMap requires non-null c argument')
        xyw = self.map(xy,c)
        return xy + (c-self.reference) * (xyw-xy)

class Compound(PixelMap):
    '''Pixel map defined as composition of other PixelMaps
    '''
    @staticmethod
    def type():
        return 'Compound'
    def __init__(self,name,elements):
        '''Create a compound PixelMap from a list of PixelMap elements
        '''
        super(PixelMap,self).__init__(name)
        self.elements = elements
        return
    def __call__(self, xy, c=None):
        # call all elements in succession
        xyw = elements[0](xy,c)
        for el in elements[1:]:
            xyw = el(xyw,c)
        return xyw

class WCS(PixelMap):
    '''Pixel map that is augmented by a projection from the world coordinates
    to definitive sky position, and toSky() will return astropy SkyCoord vectors
    corresponding to the inputs.
    Also can alternatively supply a second projection such that __call__ method 
    returns the coordinates in the second projection's system.
    A projection is anything that has toSky() and toPix() methods going to/from
    SkyCoord arrays given numpy array.
    Scale factor is multiplied into world coordinates to turn them into degrees.
    '''
    def __init__(self, name, map, projection, scale=1.):
        super(WCS,self).__init__(name)
        self.map = map
        self.projection = projection
        self.scale = scale
        self.dest = None
        return
    def reprojectTo(self, projection):
        self.dest = projection
        return
    def toSky(self, xy, c=None):
        ''' Return sky coordinates corresponding to array
        '''
        return self.projection.toSky(self.map(xy,c) * self.scale)
    def toPix(self, coords, c=None, guess=np.array([0.,0.])):
        ''' Return pixel coordinates corresponding to input SkyCoords.
        guess is a starting guess for solver, which can be either a single
        coordinate pair or an array matching coords length.
        '''
        xyw = self.projection.toPix(coords) / self.scale
        xyp = np.zeros_like(xyw) + guess  # Broadcasts to same dimensions as input
        return map.inverse(xyw, xyp, c)  # ??? tolerance???
    def __call__(self,xy,c=None):
        '''Map the input coordinates to the coordinates either in the original
        projection, or in the reprojected system if one has been given.
        '''
        xyw = self.map(xy,c)
        if self.dest is not None:
            # Execute reprojection
            coords = self.projection.toSky(xyw*self.scale)
            xyw = self.dest.toPix(coords) / self.scale
        return xyw
    
class PixelMapCollection(object):
    '''Class that holds a library of PixelMap/WCS specifications deserialized from
    a YAML file.  One can then request any PixelMap or WCS by name and be given a
    functional realization of map with that name.  Realizations are cached so that they
    are not remade every time.
    '''
    def __init__(self, filename):
        '''Create PixelMapCollection from the named YAML file
        '''
        self.root = yaml.load(open(filename))
        # Extract the WCS specifications into their own dict
        if 'WCS' in self.root:
            self.wcs = self.root.pop('WCS')
        else:
            self.wcs = {}  # Empty dictionary
        self.realizedMap = {}
        self.realizedWCS = {}
        return

    # Build a static dictionary of atomic PixelMap types
    atoms = {t.type():t for t in (Identity, Constant, Linear,
                                  Polynomial, Template, Piecewise)}

    def hasMap(self,name):
        return name in root
    def hasWCS(self,name):
        return name in wcs

    def parseAtom(self, name, **kwargs):
        '''Build a PixelMap realization for one of the atomic types
        from the specs in the kwargs dictionary
        '''
        if kwargs['Type'] not in atoms:
            raise Exception('Unknown PixelMap atomic type ' + kwargs['Type'])
        return atoms[kwargs['Type']](**kwargs)
    
    def getMap(self,name):
        ''' Return realization of PixelMap with given name
        '''
        if name in self.realizedMap:
            return self.realizedMap[name]
        if name not in root:
            raise Exception('No PixelMap with name ' + name)

        specs = self.root[name]  # The dict with specifications of this map
        if specs['Type']==ColorTerm.type():
            # Special procedure for color terms since we'll parse its captive
            map = ColorTerm(name, specs['Type'],
                            map=self.parseAtom('captive',**specs['Function']))
        elif specs['Type']==Compound.type():
            # Another special procedure for compound maps, read all its members
            elements = [self.getMap(el) for el in specs['Elements']]
            map = Compound(name, elements)
        else:
            # Map is atomic
            map = self.parseAtom(name, **specs)
        # Add new map to those already realized
        realizedMap[name] = map
        return map            

    def getWCS(self,name):
        ''' Return realization of WCS with the given name
        '''
        if name in self.realizedWCS:
            return self.realizedWCS[name]
        if name not in wcs:
            raise Exception('No WCS with name ' + name)

        specs = self.wcs[name]  # The dict with specifications of this map
        # First get the PixelMap that it modifies
        map = self.getMap(specs['MapName'])
        # Scale to turn world coords into degrees.  Note that C++ uses radians
        # for celestial coordinates and here we are in degrees by default, so
        # adjust the scale for this
        scale = specs['Scale'] * 180. / np.pi
        # Read the projection.  We recognize only 2 kinds right now:
        pspecs = specs['Projection']
        if pspecs['Type']==ICRS.type():
            projection = ICRS()
        elif pspecs['Type']==Gnomonic.type():
            ra = float(pspecs['Orientation']['RA'])
            dec = float(pspecs['Orientation']['Dec'])
            pa = float(pspecs['Orientation']['PA']) # ??? Need to check sign ???
            projection = Gnomonic(ra,dec,rotation=pa)
        else:
            raise Exception('Do not know how to parse projection of type ',pspecs['Type'])
        # Make the WCS 
        w = WCS(name, map=map, projection=projection, scale=scale)
        # Cache it
        realizedWCS[name] = w
        return w
    
########################################################
### Projections: we have just two.
########################################################

class ICRS(object):
    ''' Class giving the (trivial) projection from ICRS to ICRS coordinates, i.e.
    "pixel" coordinates are just the ICRS RA and Dec in degrees.
    '''
    def __init__(self):
        return
    def toSky(self,xy):
        return np.array(xy)
    def toPix(self,xy):
        return np.array(xy)
    
class Gnomonic(object):
    ''' Class representing a gnomonic projection about some point on
    the sky.  Can be used to go between xi,eta coordinates and SkyCoords.
    All xy units are assumed to be in degrees as are the ra, dec, and PA of
    the projection pole.
    '''
    def __init__(self, ra, dec, rotation=0.):
        '''
        Create a Gnomonic transformation by specifying the position of the
        pole (in ICRS degrees) and rotation angle of the axes relative
        to ICRS north.
        '''
        pole = co.SkyCoords(ra, dec, unit='deg',frame='ICRS')
        self.frame = pole.skyoffset_frame(rotation=co.Angle(rotation,unit='deg'))
        return
    def toSky(self,xy):
        '''
        Convert xy coordinates in the gnomonic project (in degrees) into SkyCoords.
        '''
        # First deproject the radius
        xy2d = np.atleast_2d(xy)
        rin = np.hypot(xy2d[:,0], xy2d[:,1])*np.pi/180.
        rfactor = np.where(rin>0.00001, np.arctan(rin)/rin, 1.)[:,np.newaxis]
        return co.SkyCoords(xy2d*rfactor[:,np.newaxis],unit='deg', frame=self.frame)
    def toXY(self, coords):
        '''
        Convert SkyCoord array into xy values in the gnomonic projection, in degrees
        '''
        s = coords.transform_to(self.frame)
        xy = np.vstack(s.lon.radian, s.lat.radian).T
        rin = np.hypot(xy[:,0],xy[:,1])
        rfactor = np.where(rin>0.00001, np.tan(rin) / rin, 1.) * 180. / np.pi
        return xy * rfactor[:,np.newaxis]
