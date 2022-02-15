// Test the YAMLCollector class
#include "YAMLCollector.h"
using namespace std;
using namespace astrometry;

int main(int argc,
	 char *argv[]) {
  try {
    YAMLCollector yc("i.*@testdata/dummy.yaml", "PixelMapCollection");
    YAMLCollector::Dictionary d;
    d["BAND"] = "i";
    d["DEVICE"] = "N4";
    cerr << "First add: " << yc.addMap("i/N4",d) << endl;
    d["DEVICE"] = "N5";
    cerr << "Second add: " << yc.addMap("i/N5",d) << endl;
    d["BAND"]="g";
    cerr << "Third add: " << yc.addMap("g/N5",d) << endl;
    cout << "Result:" << endl;
    cout << yc.dump() << endl;
    yc.removeMap("i/N5");
    cout << "After removal:" << endl;
    cout << yc.dump() << endl;
  } catch (std::runtime_error& e) {
    quit(e,1);
  }
}
