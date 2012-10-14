// $Id: SCAMPMap.cpp,v 1.10 2012/02/13 16:35:09 garyb Exp $
#include "SCAMPMap.h"
#include "TMV_Sym.h"

using namespace img;
using namespace poly2d;
using namespace astrometry;

// First come two functions that read the standard elements
// of the transformations from the header.

LinearMap
ReadCD(const Header* h,
       Orientation& orient) {
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
  return LinearMap(plin);
}

void
WriteCD(const LinearMap& lm,
	const Orientation& orient,
	Header& h) {
  h.replace("CTYPE1", string("RA---TAN"));
  h.replace("CTYPE2", string("DEC--TAN"));
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

SCAMPMap::SCAMPMap(const img::Header& h,
		   const Orientation* reproject): lm(0), pm(0), rm(0),
						  tp(0) {
  lm = new LinearMap(ReadCD(&h, orientIn));
  append(lm);

  tp = new TangentPlane(orientIn.getPole(), orientIn);

  Poly2d p1 = ReadPV(&h, 1);
  Poly2d p2 = ReadPV(&h, 2);
  if (p1.nCoeffs() > 1) {
    // Obtained a valid x polynomial.
    if (p2.nCoeffs() > 1) {
      // Also have valid y polynomial.  Make the map:
      pm = new PolyMap(p1, p2);
    } else {
      // Did not find PV2's.  Install default:
      Poly2d p2Identity(1);
      DVector coeffs = p2Identity.getC();
      coeffs[p2Identity.vectorIndex(0,1)] = 1.;
      p2Identity.setC(coeffs);
      pm = new PolyMap(p1, p2Identity);
    }
  } else {
    // Did not obtain any PV1's.  If there are PV2's, install
    // identity map for x coeff:
    if (p2.nCoeffs() > 1) {
      Poly2d p1Identity(1);
      DVector coeffs = p1Identity.getC();
      coeffs[p1Identity.vectorIndex(1,0)] = 1.;
      p1Identity.setC(coeffs);
      pm = new PolyMap(p1Identity, p2);
    }
  }
  if (pm) {
    // If there is a non-trivial polynomial, add it to map sequence:
    // Set 1 arcsec step for derivatives, assuming units of degrees
    pm->setPixelStep(1./3600.);
    append(pm);
  }
  // Reproject if desired:
  if (reproject) {
    orientOut.set(reproject->getPole(), reproject->getPA());
    rm = new ReprojectionMap(*tp, 
			     TangentPlane(orientOut.getPole(),orientOut),
			     DEGREE);
    append(rm);
  } else {
    orientOut = orientIn;
  }
}

SCAMPMap::~SCAMPMap() {
  if (rm) delete rm;
  if (lm) delete lm;
  if (pm) delete pm;
  if (tp) delete tp;
}

//////////////////////////////
// Now a function that will produce an image header embodying
// a linear + polynomial map that fits the input map.
//////////////////////////////
img::Header 
astrometry::FitSCAMP(Bounds<double> b,
		     const PixelMap& pm,
		     const Orientation& pmOrient,
		     const SphericalCoords& pole,
		     double tolerance) {
  const int startOrder=3;
  const int maxOrder=5;

  img::Header h;
  // The linear map from pixels to intermediate world system:
  double crval1, crval2;
  double crpix1, crpix2;
  double cd1_1, cd1_2, cd2_2, cd2_1;

  // SCAMP system will have orientation aligned to ICRS at desired pole:
  Orientation SCAMPorient(pole);
  // Make a coordinate projection map that goes from pm's world system
  // to the tangent-plane projection that we want our SCAMPMap to use
  ReprojectionMap reproject(TangentPlane(pmOrient.getPole(), pmOrient),
			    TangentPlane(pole, SCAMPorient),
			    DEGREE);
  // First thing is to set the linear map to match the
  // linearized map at center of bounds, while having
  // CRVAL at the exposure nominal center
  pole.getLonLat(crval1, crval2);
  // Degrees are standard units for FITS WCS:
  crval1 /= DEGREE;
  if (crval1<0.) crval1 += 360.;
  crval2 /= DEGREE;
  Position<double> center = b.center();
  // Map chip center to global world coords:
  double xp = center.x;
  double yp = center.y;
  double xcen, ycen;
  pm.toWorld(xp, yp, xcen, ycen);
  // Project from input map's native projection to desired projection:
  Matrix22 dWdP = reproject.dWorlddPix(xcen, ycen) * pm.dWorlddPix(xp, yp);
  reproject.toWorld(xcen, ycen, xcen, ycen);
  cd1_1 = dWdP(0,0);
  cd1_2 = dWdP(0,1);
  cd2_1 = dWdP(1,0);
  cd2_2 = dWdP(1,1);

  // CD * (xp - pix1) = xcen
  // xp - pix1 = CDinv * xcen
  // pix1 = xp - CDinv * xcen
  double det=cd1_1*cd2_2-cd1_2*cd2_1;
  crpix1 = xp - (cd2_2*xcen - cd1_2*ycen)/det;
  crpix2 = yp - (-cd2_1*xcen + cd1_1*ycen)/det;

  string s="RA---TAN";
  h.append("CTYPE1", s);
  s = "DEC--TAN";
  h.append("CTYPE2", s);
  s = "deg     ";
  h.append("CUNIT1",s);
  h.append("CUNIT2",s);
  h.append("CRPIX1", crpix1);
  h.append("CRPIX2", crpix2);
  h.append("CRVAL1", crval1);
  h.append("CRVAL2", crval2);
  h.append("CD1_1", cd1_1);
  h.append("CD1_2", cd1_2);
  h.append("CD2_1", cd2_1);
  h.append("CD2_2", cd2_2);

  // Now construct a set of test points:
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

  // starting pixel map for the first exposure with this instrument
  for (int ix=0; ix<=nx; ix++) {
    double xpix = b.getXMin() + ix*xstep;
    for (int iy=0; iy<=ny; iy++) {
      double ypix = b.getYMin() + iy*ystep;
      double xw, yw;
      pm.toWorld(xpix, ypix, xw, yw);
      // Do the change of projections:
      reproject.toWorld(xw,yw,xw,yw);
      vxp.push_back( cd1_1*(xpix-crpix1) + cd1_2*(ypix-crpix2));
      vyp.push_back( cd2_1*(xpix-crpix1) + cd2_2*(ypix-crpix2));
      vxw.push_back(xw);
      vyw.push_back(yw);
    }
  }
  // ??? select order
  int polyOrder = 3;
  double rms=0.;
  DMatrix* xcoeffs=0;
  DMatrix* ycoeffs=0;

  for (polyOrder=startOrder; polyOrder<=maxOrder; polyOrder++) {
    // Make polynomial map
    PolyMap poly(polyOrder, tolerance);

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
      poly.toWorld(vxp[i], vyp[i], xmod, ymod, derivs);
      xmod = vxw[i] - xmod;
      ymod = vyw[i] - ymod;
      for (int j=0; j<nP; j++) {
	beta[j] += xmod*derivs(0,j) + ymod*derivs(1,j);
	for (int k=0; k<=j; k++) 
	  alpha(j,k)+=derivs(0,j)*derivs(0,k) + derivs(1,j)*derivs(1,k);
      }
    }
    beta /= alpha;
    poly.setParams(beta);
    rms = 0.;
    for (int i=0; i<vxp.size(); i++) {
      double xmod, ymod;
      poly.toWorld(vxp[i], vyp[i], xmod, ymod, derivs);
      rms += pow(vxw[i] - xmod, 2.) + pow(vyw[i] - ymod,2.);
    }
    rms = sqrt(rms/vxp.size());
    //**/cerr << " rms " << rms*DEGREE/ARCSEC << " order " << polyOrder << endl;

    xcoeffs = new DMatrix(poly.getXPoly().getM());
    ycoeffs = new DMatrix(poly.getYPoly().getM());
    ycoeffs->transposeSelf();
    if (rms<tolerance) break;
  }

  if (rms>tolerance) {
    cerr << "WARNING:  SCAMPMap RMS is " << rms*DEGREE/ARCSEC
	 << " at maximum order " << polyOrder
	 << endl;
  }
  // Continue even if did not meet tolerance:
  if (polyOrder>maxOrder) polyOrder=maxOrder;

  Assert(xcoeffs->nrows()==polyOrder+1);
  Assert(xcoeffs->ncols()==polyOrder+1);
  Assert(ycoeffs->nrows()==polyOrder+1);
  Assert(ycoeffs->ncols()==polyOrder+1);
  int pvcount=0;
  for (int ord=0; ord<=polyOrder; ord++) {
    for (int yord=0; yord<=ord; yord++) {
      int xord = ord - yord;
      double coeff = (*xcoeffs)(xord, yord);
      if (coeff != 0.) {
	ostringstream oss;
	oss << "PV1_" << pvcount;
	h.append(oss.str(), coeff, "X Distortion polynomial coefficient");
      }
      coeff = (*ycoeffs)(xord, yord);
      if (coeff != 0.) {
	ostringstream oss;
	oss << "PV2_" << pvcount;
	h.append(oss.str(), coeff, "Y Distortion polynomial coefficient");
      }
      pvcount++;
      // Skip the "missing" coefficients
      if (pvcount==3) pvcount++;
      if (pvcount==11) pvcount++;
      if (pvcount==23) pvcount++;
      if (pvcount==39) pvcount++;
    }
  }
  if (xcoeffs) delete xcoeffs;
  if (ycoeffs) delete ycoeffs;
  return h;
}

void
WritePV(const PolyMap& pm,
	Header& h,
	int ipv) {
  string prefix;
  DMatrix coeffs;
  /**/cerr << "...getting coeffs" << endl;
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
  /**/cerr << "...writing coeffs" << endl;
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

Header 
SCAMPMap::writeHeader() const {
  Header h;
  if (!lm) throw AstrometryError("SCAMPMap::writeHeader() with no linear map defied");
  WriteCD(*lm, orientIn, h);
  if (pm) WritePV(*pm, h, 1);
  if (pm) WritePV(*pm, h, 2);
  return h;
}

