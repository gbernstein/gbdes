// test the FTable stuff.
#include "FTable.h"
#include <iostream>
#include <vector>
#include "Std.h"

#include "FITS.h"

using namespace std;
using namespace FITS;
using namespace img;

template <class T>
void addColumn(FTable ft, string name, long repeat) {
  if (repeat<0 || repeat>1) ft.addColumn(vector<vector<T> >(), name, repeat);
  else ft.addColumn(vector<T>(), name, repeat);
}

/*** ??? need overrides for logical & string ***/
template <class T>
void 
GetFitsColumnData(fitsfile* fptr, FTable ft, string colName, int icol, 
		  long rowStart, long rowEnd) {
  /**/cerr << "getting data for " << colName << " rowStart " << rowStart << endl;
  long nRows = rowEnd - rowStart;
  if (nRows < 1) return;
  FITS::DataType dtype = FITS::FITSTypeOf<T>();
  int status=0;
  long repeat = ft.repeat(colName);
  if (repeat < 0) {
    // Get the vector that says how long each entry is
    long repeats[nRows];
    long offsets[nRows];
    // Note adding 1 to row numbers when FITS sees them, since CFITSIO is 1-indexed:
    fits_read_descripts(fptr, icol, rowStart+1, nRows, repeats, offsets, &status);
    long maxRepeat = 0;
    for (int i=0; i<nRows; i++)
      maxRepeat = std::max(maxRepeat, repeats[i]);
    vector<T> data;
    data.reserve(maxRepeat);
    T nullValue;
    int anyNulls;
    vector<vector<T> > vv(nRows);
    for (int i=0; i<nRows; i++) {
      long length = repeats[i];
      data.resize(length);
      fits_read_col(fptr, dtype, icol, (LONGLONG) i+rowStart+1, (LONGLONG) 1, (LONGLONG) length,
		    &nullValue, &data[0], &anyNulls, &status);
      vv[i] = data;
    }
    ft.writeCells(vv, colName, rowStart);
  } else if (repeat>1) {
    long length=repeat * nRows;
    vector<T> data(length);
    T nullValue;
    int anyNulls;
    fits_read_col(fptr, dtype, icol, (LONGLONG) rowStart+1, (LONGLONG) 1, (LONGLONG) length,
		  &nullValue, &data[0], &anyNulls, &status);
    typename vector<T>::const_iterator inptr = data.begin();
    for (int i=0; i<nRows; i++) {
      ft.writeCell(vector<T>(inptr, inptr+repeat), colName, rowStart+i);
      inptr += repeat;
    }
  } else {
    // scalar column, read them all
    /**/cerr << "scalar rows: " << (LONGLONG) nRows << " sizeof element " << sizeof(T) 
	     << " dtype " << dtype << endl;
    T nullValue;
    int anyNulls;
    vector<T> data(nRows);
    fits_read_col(fptr, dtype, icol, (LONGLONG) rowStart+1, (LONGLONG) 1, (LONGLONG) nRows,
		  &nullValue, &data[0], &anyNulls, &status);/**/
    if (status) throw_CFITSIO("read scalar");
    /**/cerr << " ...read first element " << data[0];
    ft.writeCells(data, colName, rowStart);
    /**/cerr << "...done" << endl;
  }
}

int
main(int argc,
     char *argv[])
{
  try {
    string filename = argv[1];
    int hdunumber = atoi(argv[2]);

    int status=0;
    fitsfile* fitsptr=0;
    int hdutype;
    fits_open_file(&fitsptr, filename.c_str(), READWRITE, &status);
    fits_movabs_hdu(fitsptr, hdunumber, &hdutype, &status);
    if (HDUType(hdutype) != HDUBinTable) {
      cerr << "not a binary table extension" << endl;
      exit(1);
    }
    long nrows;
    int ncols;
    fits_get_num_cols(fitsptr, &ncols, &status);
    fits_get_num_rows(fitsptr, &nrows, &status);
    cout << "Table has " << ncols << " cols x " << nrows << " rows" << endl;

    img::FTable ft(nrows);

    int datatype;
    long repeat, width;
    char colName[FITS::MAX_COLUMN_NAME_LENGTH];
    vector<string> colNames;
    vector<int> colNums;
    vector<bool> colVariable;
    vector<FITS::DataType> colType;
    long bufferRows;
    fits_get_rowsize(fitsptr, &bufferRows, &status);
    cout << "Buffering suggestion: " << bufferRows << " rows." << endl;
    do {
      int colNum;
      char totallyWild[]="*";
      fits_get_colname(fitsptr, CASEINSEN, totallyWild, colName, &colNum, &status);
      if (status==0 || status==COL_NOT_UNIQUE) {
	// add to list
	colNames.push_back(colName);
	colNums.push_back(colNum);
      }
    } while (status==COL_NOT_UNIQUE);
    if (status==COL_NOT_FOUND) status=0;
    fits_clear_errmsg();
    cout << "Read " << colNames.size() << " column names." << endl;
    for (int i=0; i<colNames.size(); i++) {
      if (fits_get_eqcoltype(fitsptr, colNums[i], &datatype, &repeat, &width, &status))
	throw_CFITSIO();
      // width will give length of strings if this is string type
      bool isVariable = false;
      if (datatype < 0) {
	isVariable = true;
	datatype = -datatype;
      }
      cout << colNums[i]
	   << " " << colNames[i]
	   << " type, repeat, width: "  << datatype
	   << " " << repeat
	   << " " << width
	   << endl;

      // Add column to our table
      colVariable.push_back(isVariable);
      DataType dt = static_cast<DataType> (datatype);
      switch (dt) {
      case FITS::Tlogical:
	addColumn<bool>(ft, colNames[i], repeat);
	break;
      case FITS::Tbyte:
	addColumn<unsigned char>(ft, colNames[i], repeat);
	break;
      case FITS::Tsbyte:
	addColumn<signed char>(ft, colNames[i], repeat);
	break;
      case FITS::Tshort:
	addColumn<short>(ft, colNames[i], repeat);
	break;
      case FITS::Tushort:
	addColumn<unsigned short>(ft, colNames[i], repeat);
	break;
      case FITS::Tint32bit:
	if (FITS::CIntIsFITSLong)
	  addColumn<int>(ft, colNames[i], repeat);
	else if (!FITS::CLongIsFITSLongLong)
	  addColumn<long>(ft, colNames[i], repeat);
	else
	  throw FITSError("No intrinsic type holds 32-bit integers!");
	break;
	// ??? unsigned int possible here??? 
      case FITS::Tlonglong:
	addColumn<LONGLONG>(ft, colNames[i], repeat);
	break;
      case FITS::Tfloat:
	addColumn<float>(ft, colNames[i], repeat);
	break;
      case FITS::Tdouble:
	addColumn<double>(ft, colNames[i], repeat);
	break;
      case FITS::Tcomplex:
	addColumn<complex<float> >(ft, colNames[i], repeat);
	break;
      case FITS::Tdblcomplex:
	addColumn<complex<double> >(ft, colNames[i], repeat);
	break;
      case FITS::Tstring:
	if (isVariable) ft.addColumn(vector<string>(), colNames[i], -1);
	else {
	  Assert(repeat>=1);
	  ft.addColumn(vector<string>(), colNames[i], repeat);
	}
	// ???? Need to set up array-of-strings
	break;
      default:
	FormatAndThrow<FITSError>() << "Cannot convert data type " 
				    << dt << " to FTable Column";
      }
    }

    // Now read in the data from FITS file into our columns, using the
    // recommended row increments for buffering
    for (long rowStart = 0; rowStart < nrows; rowStart+=bufferRows) {
      long rowEnd = std::min(nrows, rowStart+bufferRows);
      for (int i=0; i<colNums.size(); i++) {
	string name = colNames[i];
	int num = colNums[i];
	/**/cerr << "Read col " << num << " <" << name << ">..." << endl;
	switch (ft.elementType(name)) {
	case FITS::Tstring:
	case FITS::Tlogical:
	  /** ??? fix later: **/ continue;
	case FITS::Tbyte:
	  GetFitsColumnData<unsigned char>(fitsptr, ft, name, num, rowStart, rowEnd);
	  break;
	case FITS::Tsbyte:
	  GetFitsColumnData<signed char>(fitsptr, ft, name, num, rowStart, rowEnd);
	  break;
	case FITS::Tshort:
	  GetFitsColumnData<short>(fitsptr, ft, name, num, rowStart, rowEnd);
	  break;
	case FITS::Tushort:
	  GetFitsColumnData<unsigned short>(fitsptr, ft, name, num, rowStart, rowEnd);
	  break;
	case FITS::Tint:
	  GetFitsColumnData<int>(fitsptr, ft, name, num, rowStart, rowEnd);
	  break;
	case FITS::Tuint:
	  GetFitsColumnData<unsigned int>(fitsptr, ft, name, num, rowStart, rowEnd);
	  break;
	case FITS::Tlong:
	  GetFitsColumnData<long>(fitsptr, ft, name, num, rowStart, rowEnd);
	  break;
	case FITS::Tulong:
	  GetFitsColumnData<unsigned long>(fitsptr, ft, name, num, rowStart, rowEnd);
	  break;
	case FITS::Tlonglong:
	  GetFitsColumnData<LONGLONG>(fitsptr, ft, name, num, rowStart, rowEnd);
	  break;
	case FITS::Tfloat:
	  GetFitsColumnData<float>(fitsptr, ft, name, num, rowStart, rowEnd);
	  break;
	case FITS::Tdouble:
	  GetFitsColumnData<double>(fitsptr, ft, name, num, rowStart, rowEnd);
	  break;
	case FITS::Tcomplex:
	  GetFitsColumnData<complex<float> >(fitsptr, ft, name, num, rowStart, rowEnd);
	  break;
	case FITS::Tdblcomplex:
	  GetFitsColumnData<complex<double> >(fitsptr, ft, name, num, rowStart, rowEnd);
	  break;
	default:
	  FormatAndThrow<FITSError>() << "Cannot convert data type " 
				      << ft.elementType(name) << " to FTable Column";
	}  // end switch on data type

      } // End column loop

    } // End row loop
    
    float f;
    string col;
    long row;
    while (cin >> col >> row) {
      ft.readCell(f, col, row);
      cout << "Value: " << f << endl;
    }

  } catch (std::runtime_error &m) {
    quit(m,1);
  }
}

