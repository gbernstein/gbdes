// Invert a big matrix
#include "UseTMV.h"
#include "Random.h"

using namespace tmv;
using namespace ran;

const int N = 20000;

int main() {
  try {
  GaussianDeviate u;
  DMatrix m(N,N);
  DVector v(N);
  for (int i=0; i<N; i++) {
      v[i] = u;
      for (int j=0; j<N; j++)
	m(i,j) = u;
  }
  cout << "Start solution" << endl;
  //  m.divideUsing(tmv::CH);
  v /= m;
  cout << "m(0,0): " << m(0,0) << " v[0] " << v[0] << endl;
  } catch (std::runtime_error& e) {
    quit(e,1);
  }
  exit(0);
}

