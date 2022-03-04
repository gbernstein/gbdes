// Random-number classes: using my old interface but replacing RNGs with
// C++ standard library random routines.
// My UniformDeviate class will also contain a bit generator, so this class
// is used where one needs either.
#ifndef RANDOM_H
#define RANDOM_H
#include <random>

//#include "Std.h"

namespace ran {
  template<class RealType=double>
  class UniformDeviate {
    // Combines a bit generator with a uniform deviate distribution.
  public:
    UniformDeviate() {seedRandom();} /*seed with value from hardware entropy source*/
    UniformDeviate(const long lseed) {seed(lseed);} /*seed with specific number*/
    RealType operator() () {return ud(gen);}
    operator RealType() {return ud(gen);}
    void seed() {seedRandom();}
    void seed(const long lseed) {gen.seed(lseed);}
    std::mt19937& generator() {return gen;}  // Expose the bit generator
  private:
    std::mt19937 gen; // The bit generator
    std::uniform_real_distribution<RealType> ud;  // The distribution
    void seedRandom() {
      std::random_device r;
      gen.seed(r());
    }
  };
  
  template<class RealType=double>
  class GaussianDeviate {
    // A unit-Gaussian deviate.  Needs a UniformDeviate underneath, either
    // making its own or sharing the one provided.
  public:
    //seed with hardware entropy
    GaussianDeviate(): ownedU(new UniformDeviate<RealType>), u(*ownedU) {}; 
    //seed with specific number:
    GaussianDeviate(const long lseed): ownedU(new UniformDeviate<RealType>(lseed)), u(*ownedU) {} 
    // Use a supplied uniform deviate:
    GaussianDeviate(UniformDeviate<RealType>& u_): ownedU(0), u(u_) {};
    ~GaussianDeviate() {if (ownedU) delete ownedU; ownedU=0;}

    RealType operator() () {return dist(u.generator());}
    operator RealType() {return dist(u.generator());}
    void seed() {u.seed();}
    void seed(const long lseed) {u.seed(lseed);}
  private:
    UniformDeviate<RealType>* ownedU;
    UniformDeviate<RealType>& u;  // Holds the bit generator
    std::normal_distribution<RealType> dist;  // The distribution
  };

  template<class RealType=double>
  class ExponentialDeviate {
    // A unit-mean exponential deviate
  public:
    //seed with hardware entropy
    ExponentialDeviate(): ownedU(new UniformDeviate<RealType>), u(*ownedU) {}; 
    //seed with specific number:
    ExponentialDeviate(const long lseed): ownedU(new UniformDeviate<RealType>(lseed)), u(*ownedU) {} 
    // Use a supplied uniform deviate:
    ExponentialDeviate(UniformDeviate<RealType>& u_): ownedU(0), u(u_) {};
    ~ExponentialDeviate() {if (ownedU) delete ownedU; ownedU=0;}

    RealType operator() () {return dist(u.generator());}
    operator RealType() {return dist(u.generator());}
    void seed() {u.seed();}
    void seed(const long lseed) {u.seed(lseed);}
  private:
    UniformDeviate<RealType>* ownedU;
    UniformDeviate<RealType>& u;  // Holds the bit generator
    std::exponential_distribution<RealType> dist;  // The distribution
  };

  template<class RealType=double>
  class BinomialDeviate {
    // A Binomial deviate for N trials each of probability p
    // The type is that of the underlying UniformDistribution, which is irrelevant
  public:
    //seed with time:
    BinomialDeviate(const int N, const double p): ownedU(new UniformDeviate<RealType>),
						    u(*ownedU),
						    dist(N,p) {}
    BinomialDeviate(const int N, const double p,
		    const long lseed): ownedU(new UniformDeviate<RealType>(lseed)),
				       u(*ownedU),
				       dist(N,p) {}
    //Use supplied uniform deviate:
    BinomialDeviate(const int N, const double p,
		    UniformDeviate<RealType>& u_): ownedU(0), u(u_),
						   dist(N,p) {}
    ~BinomialDeviate() {if (ownedU) {delete ownedU; ownedU=0;}}

    //Reset the parameters
    void reset(const int N, const double p) {dist.param(std::binomial_distribution<int>::param_type(N,p));}
    int operator() () {return dist(u.generator());}
    operator int() {return dist(u.generator());}
    void seed() {u.seed();}
    void seed(const long lseed) {u.seed(lseed);}
  private:
    UniformDeviate<RealType>* ownedU;
    UniformDeviate<RealType>& u;  // Holds the bit generator
    std::binomial_distribution<int> dist;
  };

  template<class RealType=double>
  class PoissonDeviate {
    // A Poisson deviate for process with given mean value
    // The type is that of the underlying UniformDistribution, which is irrelevant
  public:
    //seed with time:
    PoissonDeviate(const double mean): ownedU(new UniformDeviate<RealType>),
				       u(*ownedU),
				       dist(mean) {}
    PoissonDeviate(const double mean,
		   const long lseed): ownedU(new UniformDeviate<RealType>(lseed)),
				      u(*ownedU),
				      dist(mean) {}
    //Use supplied uniform deviate:
    PoissonDeviate(const double mean,
		   UniformDeviate<RealType>& u_): ownedU(0), u(u_),
						  dist(mean) {}
    ~PoissonDeviate() {if (ownedU) {delete ownedU; ownedU=0;}}

    //Reset the parameters
    void reset(const double mean) {dist.param(std::poisson_distribution<int>::param_type(mean));}
    int operator() () {return dist(u.generator());}
    operator int() {return dist(u.generator());}
    void seed() {u.seed();}
    void seed(const long lseed) {u.seed(lseed);}
  private:
    UniformDeviate<RealType>* ownedU;
    UniformDeviate<RealType>& u;  // Holds the bit generator
    std::poisson_distribution<int> dist;
  };

}  // namespace ran

#endif
