// See if simple Image & FitsImage manipulation work
#include "FitsImage.h"
using namespace img;

class Sum {
public:
  float operator()(int x, int y) const {
    return x==10 ? 200. : x+2.*y;
  }
};

int
main(int argc, char *argv[]) {
  Bounds<int> b(-99,100,-99,100);
  Image<> I(b,-1000.);

  Sum s;
  fill_pixel(I, s);

  I *= 0.5;

  cout << "Pixel 10,10: " << I(10,10) << endl;
  cout << "Pixel 20,10: " << I(20,10) << endl;
  cout << "Pixel corner: " << I(b.getXMin(), b.getYMin()) << endl;

  I.shift(1,1);

  FitsImage<> fi("test.fits", FITS::Create + FITS::OverwriteFile, -1);
  fi.copy(I);
  //  fi.retype();

  Image<> I2 = fi.extract();
  cout << "I2(1,1): " << I2(1,1) << endl;
  cout << "I2(10,10): " << I2(10,10) << endl;
  FitsImage<>::writeToFITS("test2.fits", I2,0);
  I2.setLock();
  cout << "locked I2(1,1)" << endl;
  try {
    cout << I2(1,1) << endl;
  } catch (std::runtime_error& m) {
    quit(m,1);
  }
}
