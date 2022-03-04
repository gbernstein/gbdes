// Computing Anderson-Darling statistics

#include <vector>
#include <algorithm>
#include <cmath>
#include "Random.h"

namespace stats {

  // 1-sample Anderson-Darling test.  Pass in a vector
  // of the P(x_i).
  // The return value is actually -(AD+1)*N, where AD is
  // the definition of the integral over dP.
  // Closer to zero (less negative) is better fit.
  double AD1s(vector<double>& p) {
    double Nd=p.size();
    double ad=0.;
    std::sort(p.begin(), p.end());
    if (p.front()<=0. || p.back()>=1.) return log(0.);
    for (int i=0; i<p.size(); i++) {
      ad += (2*i+1)*log(p[i]) + (2*p.size()-2*i-1)*log(1.-p[i]);
    }
    
    return (ad/p.size());
  }

  // Calculate AD1s for one randomly generated sample of size N
  double randomAD1s(int N, ran::UniformDeviate& u) {
    vector<double> p(N);
    for (int i=0; i<N; i++)
      p[i]=u();
    return AD1s(p);
  }

  // Calculate probability of AD1s>= x for vector length N with
  // nSamp random samples
  double pAD1s(double x, int N, int nSamp=1000) {
    int nGE=0;
    ran::UniformDeviate u;
    for (int i=0; i<nSamp; i++)
      if (randomAD1s(N,u)>=x) nGE++;
    double f=nGE;
    return f/nSamp;
  }

  // K-sample Anderson-Darling
  double ADks(vector<vector<double>* > vv) {
    int nSamples = vv.size();
    
    vector<double> merged;
    for (int iSample=0; iSample<nSamples; iSample++) {
      // Sort each sample and make a merged sorted list
      std::sort(vv[iSample]->begin(), vv[iSample]->end());
      merged.reserve(merged.size() + vv[iSample]->size());
      vector<double>::iterator mm=merged.end();
      merged.insert(merged.end(), 
		    vv[iSample]->begin(), vv[iSample]->end());
      std::inplace_merge(merged.begin(),
			 mm,
			 merged.end());
      //cerr << "After merge " << iSample << ": " << endl;
      //for (int j=0; j<merged.size(); j++) cerr << merged[j] << " " ;
      //cerr << endl;
    }

    double sumAD = 0.;
    double invN = 1./merged.size();
    for (int iSample=0; iSample<nSamples; iSample++) {
      // Make p vector for each sample & get AD
      vector<double> p(*(vv[iSample]));
      int lessIndex=0;  //first merged item >= this one
      int gtrIndex=0;    //first merged > this one.
      if (p.empty()) continue;
      for (int i=0; i<p.size(); i++) {
	double x=p[i];
	while ( (merged[lessIndex] < x)
		&& lessIndex<merged.size()) lessIndex++;
	while (merged[gtrIndex] <= x) {
	  if (gtrIndex>=merged.size()) break;
	  gtrIndex++;
	}
	p[i] = invN*(lessIndex + 0.5*(gtrIndex-lessIndex));
      }
      sumAD += AD1s(p);
    }

    return sumAD;
  }

  double randomADks(vector<int>& nn, ran::UniformDeviate& u) {
    vector<vector<double>*> vv(nn.size());

    for (int isamp=0; isamp<nn.size(); isamp++) {
      // make fake data
      vector<double>* v=new vector<double>(nn[isamp]);
      for (int i=0; i<v->size(); i++)
	(*v)[i]=u();
      vv[isamp] = v;
    }
    double result = ADks(vv);

    for (int isamp=0; isamp<vv.size(); isamp++)
      delete vv[isamp];

    return result;
  }

  double pADks(double x, vector<int> nn, int nSamp=1000) {
    int nGE=0;
    ran::UniformDeviate u;
    for (int i=0; i<nSamp; i++)
      if (randomADks(nn,u)>=x) nGE++;
    double f=nGE;
    return f/nSamp;
  }

  double pADks(double x, vector<vector<double>*> vv, int nSamp=1000) {
    vector<int> nn(vv.size());
    for (int i=0; i<nn.size(); i++) nn[i] = vv[i]->size();
    return pADks(x,nn,nSamp);
  }
} //namespace stats

