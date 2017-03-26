// Code for PhotoPrior class that pushes different exposures to common zeropoint/airmass/color terms
// or just pushes these to predetermined values to break degeneracies.
#include "PhotoMatch.h"
#include <map>
#include <StringStuff.h>

using namespace photometry;

PhotoPrior::PhotoPrior(list<PhotoPriorReferencePoint> points_,
		       double sigma_,
		       string name_,
		       double zeropoint,
		       double airmassCoefficient,
		       double colorCoefficient,
		       bool freeZeropoint,
		       bool freeAirmass,
		       bool freeColor): sigma(sigma_),
					points(points_),
					name(name_),
					nFree(0),
					globalStartIndex(-1),
					globalMapNumber(-1),
					m(zeropoint),
					a(airmassCoefficient),
					b(colorCoefficient),
					mIsFree(freeZeropoint),
					aIsFree(freeAirmass),
					bIsFree(freeColor)
{
  if (mIsFree) nFree++;
  if (aIsFree) nFree++;
  if (bIsFree) nFree++;
  countFit();
}

DVector
PhotoPrior::getParams() const {
  DVector p(nFree);
  if (nFree==0) return p;
  int index = 0;
  if (mIsFree) p[index++] = m;
  if (aIsFree) p[index++] = a;
  if (bIsFree) p[index++] = b;
  return p;
}

void
PhotoPrior::setParams(const DVector& p) {
  Assert(p.size() == nFree);
  if (nFree==0) return;
  int index=0;
  if (mIsFree) m = p[index++];
  if (aIsFree) a = p[index++];
  if (bIsFree) b = p[index++];
}

void
PhotoPrior::countFit() {
  nFit = 0;
  for (list<PhotoPriorReferencePoint>::const_iterator i = points.begin();
       i != points.end();
       ++i) {
    if (!(i->isClipped)) nFit++;
  }
}

double
PhotoPrior::chisq(int& dof) const {
  // We are not going to use the prior in chisq if it would be degenerate:
  if (isDegenerate()) return 0.;

  double chi = 0.;
  for (list<PhotoPriorReferencePoint>::const_iterator i = points.begin();
       i != points.end();
       ++i) {
    if (i->isClipped) continue;
    double dm = i->magOut - m - a*(i->airmass-1.) - b*i->args.color;
    chi += dm * dm;
  }
  dof += nFit - nFree;
  return chi / (sigma*sigma);
}

void
PhotoPrior::remap() {
  for (list<PhotoPriorReferencePoint>::iterator i = points.begin();
       i != points.end();
       ++i) 
    i->magOut = i->map->forward(i->magIn, i->args);
}
  
int
PhotoPrior::accumulateChisq(double& chisq,
			    DVector& beta,
			    astrometry::AlphaUpdater& updater) {
  if (isDegenerate()) return 0; // Can't use this Prior if it's degenerate.

  int nP = beta.size();
  double wt = 1./(sigma*sigma);

  // Update mapping and save derivatives for each detection:
  for (list<PhotoPriorReferencePoint>::iterator i=points.begin();
       i!=points.end(); 
       ++i) {
    if (i->isClipped) continue;
    int npi = i->map->nParams();
    DVector derivs1(npi);
    if (npi>0) {
      i->magOut = i->map->forwardDerivs(i->magIn, i->args, derivs1);
    } else {
      i->magOut = i->map->forward(i->magIn, i->args);
    }
    double dm = i->magOut - m - a*(i->airmass-1.) - b*i->args.color;

    // Here are the derivs with respect to the prior's parameters:
    DVector derivs2(nParams());
    int index = 0;
    if (mIsFree) derivs2[index++] = -1.;
    if (aIsFree) derivs2[index++] = -(i->airmass-1.);
    if (bIsFree) derivs2[index++] = -i->args.color;
      
    chisq += dm * dm * wt;

    // Accumulate derivatives:
    int istart=0;  // Index into this point's derivative vector
    for (int iMap=0; iMap<i->map->nMaps(); iMap++) {
      int ip=i->map->startIndex(iMap);
      int np=i->map->nSubParams(iMap);
      if (np==0) continue;
      int mapNumber1 = i->map->mapNumber(iMap);

      // Derivs for just this part of the submap:
      tmv::VectorView<double> sub1=derivs1.subVector(istart,istart+np);
      beta.subVector(ip, ip+np) -= wt*dm*sub1;

      // Augment the portion of alpha straddling diagnonal:
      //**      alpha.subSymMatrix(ip, ip+np) += wt*(sub1 ^ sub1));
      updater.rankOneUpdate(mapNumber1, ip, sub1, 
			    mapNumber1, ip, sub1, wt);

      // Augment the alpha from cross-derivatives with the prior parameters 
      if (nParams() > 0) {
	int mapNumber2 = mapNumber();
	updater.rankOneUpdate(mapNumber2, startIndex(), derivs2.view(),
			      mapNumber1, ip, sub1, wt);
      }

      // Augment alpha with cross-derivatives from other maps' params:
      int istart2 = istart+np;
      for (int iMap2=iMap+1; iMap2<i->map->nMaps(); iMap2++) {
	int ip2=i->map->startIndex(iMap2);
	int np2=i->map->nSubParams(iMap2);
	int mapNumber2 = i->map->mapNumber(iMap2);
	if (np2==0) continue;
	tmv::VectorView<double> sub2=derivs1.subVector(istart2,istart2+np2);
	// Note that subMatrix here will not cross diagonal.
	//**	  alpha.subMatrix(ip, ip+np, ip2, ip2+np2) += wt*(sub1 ^ sub2);
	updater.rankOneUpdate(mapNumber2, ip2, sub2,
			      mapNumber1, ip,  sub1, wt);
	istart2+=np2;
      }
      istart+=np;
    } // outer parameter segment loop

    // Now add derivs wrt prior's parameters to beta and alpha:
    if (nParams() > 0) {
      int ip = startIndex();
      int np = nParams();
      beta.subVector(ip, ip+np) -= wt*dm*derivs2;

      // Augment the portion of alpha straddling diagnonal:
      //**      alpha.subSymMatrix(ip, ip+np) += wt*(derivs2 ^ derivs2));
      updater.rankOneUpdate(mapNumber(), ip, derivs2.view(), 
			    mapNumber(), ip, derivs2.view(), wt);
    }
  } // reference point loop

  return nFit - nFree;
}

bool
PhotoPrior::sigmaClip(double sigThresh) {
  list<PhotoPriorReferencePoint>::iterator worst=points.end();
  double maxDev = sigThresh*sigma;
  for (list<PhotoPriorReferencePoint>::iterator i=points.begin();
       i!=points.end(); 
       ++i) {
    if (i->isClipped) continue;
    double dm = i->magOut - m - a*(i->airmass-1.) - b*i->args.color;
    if (abs(dm) > maxDev) {
      worst = i;
      maxDev = abs(dm);
    }
  }

  if (worst == points.end()) return false;

  // Have one to clip:
  worst->isClipped = true;
  nFit--;
  // Will also clip other Points from the same exposure (different color)
  for (list<PhotoPriorReferencePoint>::iterator i=points.begin();
       i!=points.end(); 
       ++i) {
    if (i==worst || i->isClipped) continue;
    if (i->exposureName == worst->exposureName) {
      i->isClipped = true;
      nFit--;
  }  
  }
  return true;
}

void
PhotoPrior::clipAll() {
  for (list<PhotoPriorReferencePoint>::iterator i=points.begin();
       i!=points.end(); 
       ++i)
    i->isClipped = true;

  nFit = 0;
}

void
PhotoPrior::reportHeader(ostream& os) {
  os <<  
    "Name        chisq  / dof  sigma   zpt    airmass  color   \n"
    "* Exposure   Device   Residual  magIn  airmass  color   \n"
    "----------------------------------------------------------"
     << endl;
}

void
PhotoPrior::report(ostream& os) const {
  stringstuff::StreamSaver ss(os);

  double chi = 0.;
  int dof;

  if (isDegenerate()) {
    dof = nFit - nFree;
  } else {
    chi = chisq(dof);
  }
  os << setw(12) << left << getName()
     << " "  << fixed << right << noshowpos << setprecision(1) << setw(6) << chi
     << " / " << setw(3) << dof
     << " " << setprecision(3) << setw(5) << getSigma()
     << " " << showpos << setprecision(3) << setw(7) << getZeropoint()
     << " " << (zeropointIsFree() ? " " : "*" )
     << " " << setprecision(3) << setw(6) << getAirmass()
     << " " << (airmassIsFree() ? " " : "*" )
     << " " << setprecision(3) << setw(6) << getColor()
     << " " << (colorIsFree() ? " " : "*" )
     << endl;

  for (list<PhotoPriorReferencePoint>::const_iterator i = points.begin();
       i != points.end();
       ++i) {
    double resid = i->magOut - (m + a*(i->airmass-1.) + b*i->args.color);
    os << (i->isClipped ? "- " : "+ ")
       << setw(12) << left << i->exposureName
       << " " << setw(8) << i->deviceName
       << " " << right << showpos << setprecision(3) << setw(6) << resid
       << " " << noshowpos << setprecision(3) << setw(7) << i->magIn
       << " " << setprecision(3) << setw(5) << i->airmass
       << " " << showpos << setprecision(3) << setw(6) << i->args.color
       << endl;
  }
}
