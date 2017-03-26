// Prepare an input file for WCSFoF
#include "FitsTable.h"
#include <sstream>

using namespace std;
using namespace FITS;
using namespace img;


int
main(int argc,
     char *argv[])
{

  try {
    int nfiles = 6;
    FTable tab(nfiles);
    vector<string> bases(nfiles);
    bases[0]="2MASS_2239-0544_r44";
    bases[1]="I1110738p";
    bases[2]="I1110739p";
    bases[3]="I1110740p";
    bases[4]="I1110741p";
    bases[5]="I1110742p";

    tab.addColumn(bases, "Exposure");
    string dummy="2MASS";
    tab.writeCell(dummy, "Exposure", 0);

    vector<string> names=bases;
    for (int i=0; i<nfiles; i++) 
      names[i] += ".cat";
    tab.addColumn(names, "Filename");
    
    names = bases;
    for (int i=0; i<nfiles; i++) 
      names[i] += ".head";
    names[0] = "ICRS";
    tab.addColumn(names, "WCSFILE");
    

    vector<string> fields(nfiles,"@_NEAREST");
    fields[0]="Wegner";
    tab.addColumn(fields, "Field");

    vector<string> ra(nfiles,"@RA");
    ra[0] = "0";
    /*
    ra[1] = "22:39:33.00";
    ra[2] = "22:39:32.00";
    ra[3] = "22:39:33.30";
    ra[4] = "22:39:34.00";
    ra[5] = "22:39:32.70";  */

    vector<string> dec(nfiles,"@DEC");
    dec[0]="0";
    /*
    dec[1] = "-5:44:56.0";
    dec[2] = "-5:45:00.5";
    dec[3] = "-5:45:11.0";
    dec[4] = "-5:44:51.5";
    dec[5] = "-5:44:41.0";*/
    tab.addColumn(ra,"RA");
    tab.addColumn(dec,"DEC");

    vector<string> instr(nfiles,"CFHTi");
    instr[0] = "Reference";
    tab.addColumn(instr,"Instrument");

    vector<string> affinity(nfiles,"I");
    tab.addColumn(affinity,"Affinity");

    vector<string> device(nfiles,"@CCDNICK");
    tab.addColumn(device,"DEVICE");

    vector<string> xkey(nfiles,"XWIN_IMAGE");
    xkey[0] = "X_WORLD";
    tab.addColumn(xkey,"XKey");

    vector<string> ykey(nfiles,"YWIN_IMAGE");
    ykey[0] = "Y_WORLD";
    tab.addColumn(ykey,"YKey");

    vector<string> idkey(nfiles,"NUMBER");
    idkey[0] = "@_ROW";
    tab.addColumn(idkey,"IDKey");

    vector<string> select(nfiles,"FLAGS==0 && ERRAWIN_IMAGE<0.1");
    select[0] = "ERRA_WORLD < 0.5/3600.";
    tab.addColumn(select,"Select");

    vector<int> ext(nfiles,-1);
    tab.addColumn(ext,"Extension");

  FitsTable file("index.fits",FITS::ReadWrite + FITS::Create + FITS::OverwriteFile);
  file.copy(tab);
  } catch (std::runtime_error &m) {
    quit(m,1);
  }
}
