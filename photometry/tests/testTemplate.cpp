// See if the template photo maps work as expected

#include <iostream>
#include "PhotoTemplate.h"
#include "Std.h"
#include "yaml-cpp/yaml.h"
using namespace photometry;
using namespace std;

int
main(int argc, char *argv[])
{
  try {
    PhotoTemplate test(argv[1], argv[2], "TemplateTest");
    test.setParams(DVector(1, 0.5));
    cout << "parameter: " << test.getParams() << endl;
    PhotoArguments args;
    while (cin >> args.xDevice >> args.yDevice) {
      cout << "Result: " << test.forward(0., args) << endl;
      DVector derivs(1);
      test.forwardDerivs(0., args, derivs);
      cout << "Deriv: " << derivs[0] << endl;
      cout << "Next: ";
    }
    YAML::Emitter out;
    test.write(out);
    cout << out.c_str() << endl;
  } catch (std::runtime_error& m) {
    quit(m,1);
  }
}


