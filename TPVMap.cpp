#include "TPVMap.h"
#include "TMV_Sym.h"
#include "PixelMapCollection.h"

using namespace img;
using namespace poly2d;
using namespace astrometry;

// First come two functions that read the standard elements
// of the transformations from the header.

LinearMap
ReadCD(const Header* h,
       Orientation& orient,
       string name) {
  double lon, lat;
  double crpix1, crpix2;
  double cd11, cd12, cd21, cd22;
  string ctype1,ctype2;

  DVector plin2((size_t) 6, 0.);
  if (!h->getValue("CTYPE1", ctype1))
    throw AstrometryError("ReadCD could not find CTYPE1");
  if (!h->getValue("CTYPE2", ctype2))
    throw AstrometryError("ReadCD could not find CTYPE2");
  if (!h->getValue("CRPIX1", crpix1))
    throw AstrometryError("ReadCD could not find CRPIX1");
  if (!h->getValue("CRPIX2", crpix2))
    throw AstrometryError("ReadCD could not find CRPIX2");
  if (!h->getValue("CRVAL1", lon))
    throw AstrometryError("ReadCD could not find CRVAL1");
  if (!h->getValue("CRVAL2", lat))
    throw AstrometryError("ReadCD could not find CRVAL2");
  if (!h->getValue("CD1_1", cd11))
    throw AstrometryError("ReadCD could not find CD1_1");
  if (!h->getValue("CD1_2", cd12))
    throw AstrometryError("ReadCD could not find CD1_2");
  if (!h->getValue("CD2_1", cd21))
    throw AstrometryError("ReadCD could not find CD2_1");
  if (!h->getValue("CD2_2", cd22))
    throw AstrometryError("ReadCD could not find CD2_2");

  if ( (ctype1!="RA---TAN" || ctype2!="DEC--TAN")
       && (ctype1!="RA---TPV" || ctype2!="DEC--TPV")
       && (ctype1!="RA---TNX" || ctype2!="DEC--TNX"))
    FormatAndThrow<AstrometryError>() << "ReadCD does not know CTYPEs " 
				      << ctype1 << " or " << ctype2;
  // *** ??? Note assumption that CRVAL is ICRS coordinate in degrees
  SphericalICRS pole(lon*DEGREE, lat*DEGREE);
  orient.set(pole,0.);

  // Get parameters of the linear map
  DVector plin(6);
  plin[0] = -cd11*crpix1 - cd12*crpix2;
  plin[1] = cd11;
  plin[2] = cd12;
  plin[3] = -cd21*crpix1 - cd22*crpix2;
  plin[4] = cd21;
  plin[5] = cd22;
  return LinearMap(plin, name);
}

void
WriteCD(const LinearMap& lm,
	const Orientation& orient,
	Header& h) {
  h.replace("CTYPE1", string("RA---TPV"));
  h.replace("CTYPE2", string("DEC--TPV"));
  double ra, dec;
  SphericalICRS pole = orient.getPole();
  pole.getLonLat(ra,dec);
  ra /= DEGREE;
  if (ra < 0.) ra += 360.;
  dec /= DEGREE;
  h.replace("CRVAL1", ra);
  h.replace("CRVAL2", dec);

  double xpix, ypix;
  lm.toPix(0., 0., xpix, ypix);
  h.replace("CRPIX1", xpix);
  h.replace("CRPIX2", ypix);
  Matrix22 m = lm.dWorlddPix(xpix, ypix);
  h.replace("CD1_1", m(0,0));
  h.replace("CD1_2", m(0,1));
  h.replace("CD2_1", m(1,0));
  h.replace("CD2_2", m(1,1));
}

Poly2d
ReadPV(const Header* h, int ipv) {
  string prefix;
  if (ipv==1) 
    prefix="PV1_";
  else if (ipv==2) 
    prefix="PV2_";
  else 
    throw AstrometryError("SCAMPMap::ReadPV needs ipv=1 or 2");
  vector<double> coeffs;
  // Convention is that PV1_1 and PV2_1 default to 1. while other
  // coefficients default to 0, so that absence of any keywords equals
  // identity transformation. Keep track of whether linear term PVx_1 was
  // given:
  bool foundLinear = false;
  h->rewind();
  for (h->rewind(); !h->atEnd(); h->incr()) {
    string key = h->current()->getKeyword();
    if (key.compare(0, prefix.size(), prefix)!=0) continue;
    key.erase(0,prefix.size());
    int index = atoi(key.c_str());
    if (index<0) 
      FormatAndThrow<AstrometryError>() 
	<< "Remainder of PV[12]_x key is not integer: " << key;
    if (index==1) foundLinear = true;
    const HdrRecord<double>* p = dynamic_cast<const HdrRecord<double>*> (h->current());
    if (!p) 
      throw AstrometryError("Could not cast PV1_* record to double");

    double v= p->Value();
    // Have index & value, add to a vector
    // Magic: 3, 11, 23, and 39 are missing, they are (non-analytic) odd powers of r
    // in some documentation.
    if (index==3 || index==11 || index==23 || index==39) continue;
    if (index>=39) index--;
    if (index>=23) index--;
    if (index>=11) index--;
    if (index>=3) index--;
    while (index>=coeffs.size()) coeffs.push_back(0.);
    Assert(index<coeffs.size());
    coeffs[index]=v;
  }

  // If no keywords were found, return a zero-order polynomial
  if (coeffs.empty()) return Poly2d(0);

  // Now convert coefficient vector into a matrix
  // What is the order of polynomial? Make it linear, at least
  int order=1;
  while ( (order+1)*(order+2)/2 < coeffs.size()) order++;
  DMatrix cm(order+1, order+1, 0.);
  int i=0; int j=0;
  for (int k=0; k<coeffs.size(); k++) {
    Assert(i>=0 && j>=0 && i<=order && j<=order);
    cm(i,j) = coeffs[k];
    if (i==0) {
      i=j+1; j=0;
    } else {
      i--; j++;
    }
  }
  // Put in default value for linear term if none was given:
  if (!foundLinear) cm(1,0) = 1.;

  // Transpose for PV2:
  if (ipv==2) cm.transposeSelf();

  // Then make a Poly2d from it
  return Poly2d(cm);
}

void
WritePV(const PolyMap& pm,
	 Header& h,
	int ipv) {
  string prefix;
  DMatrix coeffs;
  if (ipv==1) {
    prefix="PV1_";
    DMatrix cm = pm.getXPoly().getM();
    coeffs.resize(cm.nrows(), cm.ncols());
    coeffs = cm;
  } else if (ipv==2) {
    prefix="PV2_";
    DMatrix cm = pm.getYPoly().getM();
    coeffs.resize(cm.nrows(), cm.ncols());
    coeffs = cm;
    coeffs.transposeSelf();
  } else {
    throw AstrometryError("SCAMPMap::WritePV needs ipv=1 or 2");
  }
  int order = coeffs.nrows()-1;
  int i=0;
  int j=0;
  int pvnumber=0;
  while (i+j <= order) {
    if (coeffs(i,j) != 0.) {
      ostringstream oss;
      oss << prefix << pvnumber;
      h.replace(oss.str(), coeffs(i,j));
    }
    // Move to next PV number, skipping the nonsense ones:
    pvnumber++;
    if (pvnumber==3) pvnumber++;
    if (pvnumber==11) pvnumber++;
    if (pvnumber==23) pvnumber++;
    if (pvnumber==39) pvnumber++;
    // Move down diagonals of matrix:
    if (i==0) {
      i=j+1; j=0;
    } else {
      i--; j++;
    }
  }
}

/////////////////////////////////////////
// Now the higher-level read/write from headers to Wcs's.
/////////////////////////////////////////

// step for derivatives of polynomials is 1 arcsec, in degrees
const double POLYSTEP = 1./3600.;
// Suffixes to give for the linear and polynomial maps' names:
const string linearSuffix = "_cd";
const string polySuffix = "_pv";

Wcs* 
astrometry::readTPV(const img::Header& h, string name) {
  Orientation orientIn;
  LinearMap lm = ReadCD(&h, orientIn, name+linearSuffix);

  Gnomonic tp(orientIn.getPole(), orientIn);
  PixelMap* pv=0;
  string polyName = name + polySuffix;

  Poly2d p1 = ReadPV(&h, 1);
  Poly2d p2 = ReadPV(&h, 2);
  if (p1.nCoeffs() > 1) {
    // Obtained a valid x polynomial.
    if (p2.nCoeffs() > 1) {
      // Also have valid y polynomial.  Make the map:
      pv = new PolyMap(p1, p2, polyName,POLYSTEP);
    } else {
      // Did not find PV2's.  Install default:
      Poly2d p2Identity(1);
      DVector coeffs = p2Identity.getC();
      coeffs[p2Identity.vectorIndex(0,1)] = 1.;
      p2Identity.setC(coeffs);
      pv = new PolyMap(p1, p2Identity,polyName,POLYSTEP);
    }
  } else {
    // Did not obtain any PV1's.  If there are PV2's, install
    // identity map for x coeff:
    if (p2.nCoeffs() > 1) {
      Poly2d p1Identity(1);
      DVector coeffs = p1Identity.getC();
      coeffs[p1Identity.vectorIndex(1,0)] = 1.;
      p1Identity.setC(coeffs);
      pv = new PolyMap(p1Identity, p2, polyName, POLYSTEP);
    }
  }
  list<PixelMap*> pmlist;
  pmlist.push_back(&lm);
  if (pv) pmlist.push_back(pv);
  // Create a SubMap that owns a duplicate of these one or two maps (no sharing):
  SubMap sm(pmlist, name, false);
  // Delete the polymap if we made it:
  if (pv) delete pv;
  // Return the Wcs, which again makes its own copy of the maps (no sharing):
  return new Wcs(&sm, tp, name, DEGREE, false);
}

Header
astrometry::writeTPV(const Wcs& w) {
  Header h;
  const LinearMap* lm=0;
  const PolyMap* pv=0;
  const PixelMap* pm = w.getMap();
  if ( dynamic_cast<const LinearMap*> (pm)) {
    // The PixelMap is purely linear:
    lm =  dynamic_cast<const LinearMap*> (pm);
  } else {
    const SubMap* sm = dynamic_cast<const SubMap*> (pm);
    if (!sm)
      throw AstrometryError("writeTPV(): input Wcs <" + w.getName()
			    + "> has wrong type of PixelMap <" + sm->getName() + ">");
    if (sm->nMaps() < 1 || sm->nMaps() > 2)
      throw AstrometryError("writeTPV(): input Wcs <" + w.getName()
			    + "> has SubMap <" + sm->getName()
			    + "> with wrong number of PixelMap components");
    lm = dynamic_cast<const LinearMap*> (sm->getMap(0));
    if (!lm)
      throw AstrometryError("writeTPV(): First component of input Wcs <" + w.getName()
			    + "> is not LinearMap");
    if (sm->nMaps()>1) {
      pv = dynamic_cast<const PolyMap*> (sm->getMap(1));
      if (!pv)
	throw AstrometryError("writeTPV(): Second component of input Wcs <" + w.getName()
			      + "> is not PolyMap");
    }
  }

  const Gnomonic* gn = dynamic_cast<const Gnomonic*> (w.getNativeCoords());
  if (!gn)
    throw AstrometryError("writeTPV(): input Wcs <" + w.getName() 
			  + "> does not have Gnomonic projection.");
  if (gn->getOrient()->getPA() != 0.)
    throw AstrometryError("writeTPV() Orientation with non-zero PA not implemented");

  WriteCD(*lm, *gn->getOrient(), h);
  if (pv) WritePV(*pv, h, 1);
  if (pv) WritePV(*pv, h, 2);
  return h;
}

//////////////////////////////
// Now a function that will fit a linear + polynomial + gnomonic deprojection
// to an arbitrary Wcs
//////////////////////////////
Wcs*
astrometry::fitTPV(Bounds<double> b,
		   const Wcs& wcsIn,
		   const SphericalCoords& tpvPole,
		   string name,
		   double tolerance) {
  const int startOrder=3;
  const int maxOrder=5;

  const string linearName = name + linearSuffix;
  const string polyName = name + polySuffix;

  img::Header h;
  // The linear map from pixels to intermediate world system:
  double crval1, crval2;
  double crpix1, crpix2;
  double cd1_1, cd1_2, cd2_2, cd2_1;

  // Coordinates in the desired TPV gnomonic projection:
  Orientation tpvOrient(tpvPole);
  Gnomonic tpvCoords(tpvOrient);

  // Make the LinearMap that has the same behavior near the center of bounds:
  //  xworld = v[0] + v[1]*xpix + v[2]*ypix;
  //  yworld = v[3] + v[4]*xpix + v[5]*ypix;
  DVector v(6);

  Position<double> center = b.center();
  // Map chip center and deviates in x and y to coords in the desired gnomonic projection:
  double xp = center.x;
  double yp = center.y;
  tpvCoords.convertFrom(wcsIn.toSky(xp, yp));
  double x0, y0;
  tpvCoords.getLonLat(x0,y0);
  x0 /= DEGREE;  y0 /= DEGREE;
  xp = center.x + wcsIn.getPixelStep();
  tpvCoords.convertFrom(wcsIn.toSky(xp, yp));
  double dx, dy;
  tpvCoords.getLonLat(dx,dy);
  v[1] = (dx/DEGREE - x0)/wcsIn.getPixelStep();
  v[4] = (dy/DEGREE - y0)/wcsIn.getPixelStep();

  xp = center.x;
  yp = center.y + wcsIn.getPixelStep();
  tpvCoords.convertFrom(wcsIn.toSky(xp, yp));
  tpvCoords.getLonLat(dx,dy);
  v[2] = (dx/DEGREE - x0)/wcsIn.getPixelStep();
  v[5] = (dy/DEGREE - y0)/wcsIn.getPixelStep();

  v[0] = x0 - (v[1]*center.x + v[2]*center.y);
  v[3] = y0 - (v[4]*center.x + v[5]*center.y);

  LinearMap lm(v, linearName);

  // Now construct a set of test points to which we'll fit polynomial:
  const int nGridPoints=400;	// Number of test points for map initialization
  vector<double> vxp;
  vector<double> vyp;
  vector<double> vxw;
  vector<double> vyw;

  double step = sqrt(b.area()/nGridPoints);
  int nx = static_cast<int> (ceil((b.getXMax()-b.getXMin())/step));
  int ny = static_cast<int> (ceil((b.getYMax()-b.getYMin())/step));
  double xstep = (b.getXMax()-b.getXMin())/nx;
  double ystep = (b.getYMax()-b.getYMin())/ny;

  // Find post-linear coords of each test point and desired sky posn in target coords
  for (int ix=0; ix<=nx; ix++) {
    for (int iy=0; iy<=ny; iy++) {
      double xpix = b.getXMin() + ix*xstep;
      double ypix = b.getYMin() + iy*ystep;
      double xw, yw;
      tpvCoords.convertFrom(wcsIn.toSky(xpix, ypix));
      tpvCoords.getLonLat(xw, yw);
      xw /= DEGREE;
      yw /= DEGREE;
      lm.toWorld(xpix, ypix, xpix, ypix);
      vxp.push_back(xpix);
      vyp.push_back(ypix);
      vxw.push_back(xw);
      vyw.push_back(yw);
    }
  }

  // select order
  int polyOrder = 3;
  double rms=0.;

  PolyMap* pv=0;
  for (polyOrder=startOrder; polyOrder<=maxOrder; polyOrder++) {
    // Make polynomial map
    PolyMap poly(polyOrder, polyName, tolerance);

    // Fit map to test points
    int nP = poly.nParams();
    DVector beta(nP, 0.);
    tmv::SymMatrix<double> alpha(nP, 0.);
    DMatrix derivs(2,nP);
    Assert(vxp.size()==vyp.size());
    Assert(vxp.size()==vxw.size());
    Assert(vxp.size()==vyw.size());
    for (int i=0; i<vxp.size(); i++) {
      double xmod, ymod;
      poly.toWorldDerivs(vxp[i], vyp[i], xmod, ymod, derivs);
      for (int j=0; j<nP; j++) {
	beta[j] += vxw[i]*derivs(0,j) + vyw[i]*derivs(1,j);
	for (int k=0; k<=j; k++) 
	  alpha(j,k)+=derivs(0,j)*derivs(0,k) + derivs(1,j)*derivs(1,k);
      }
    }
    beta /= alpha;
    poly.setParams(beta);
    if (pv) delete pv;
    pv = new PolyMap(poly);
    rms = 0.;
    for (int i=0; i<vxp.size(); i++) {
      double xmod, ymod;
      poly.toWorld(vxp[i], vyp[i], xmod, ymod);
      rms += pow(vxw[i] - xmod, 2.) + pow(vyw[i] - ymod,2.);
    }
    rms = sqrt(rms/vxp.size());
    //**/cerr << " rms " << rms*DEGREE/ARCSEC << " order " << polyOrder << endl;
    if (rms<tolerance) break;
  }

  if (rms>tolerance) {
    cerr << "WARNING:  SCAMPMap RMS is " << rms*DEGREE/ARCSEC
	 << " at maximum order " << polyOrder-1
	 << " fitting Wcs " << wcsIn.getName()
	 << endl;
  }

  // Continue even if did not meet tolerance:
  list<PixelMap*> pmlist;
  pmlist.push_back(&lm);
  pmlist.push_back(pv);
  // Compound the two maps into a SubMap, do not share 
  SubMap sm(pmlist, name, false);
  delete pv;
  // And return a Wcs built from this SubMap (Wcs owns the maps)
  return new Wcs(&sm, tpvCoords, name, DEGREE, false);
}



