// Read/Write FTables from FITS extensions

#include "FITSTable.h"

using namespace img;
using namespace FITS;

FITSTable::FITSTable(string filename, int hduNumber_): fptr(0), 
						       fname(filename),
						       hduNumber(hduNumber_) {
  int status=0;
  int hdutype;
  fits_open_file(&fptr, fname.c_str(), READONLY, &status);
  fits_movabs_hdu(fptr, hduNumber, &hdutype, &status);
  if ( (HDUType(hdutype) != HDUBinTable) ||
       (HDUType(hdutype) != HDUAsciiTable) )
    throw FITSError("Not a FITS table extension");

  fits_get_num_cols(fptr, &ncols, &status);
  fits_get_num_rows(fptr, &nrows, &status);
  if (status) throw_CFITSIO("Opening FITSTable in file " + fname);
}

void
FITSTable::moveTo() const {
  int status=0;
  int hdutype;
  fits_movabs_hdu(fptr, hduNumber, &hdutype, &status);
  if (status) throw_CFITSIO();
}

~FITSTable() {
  if (!fptr) return;
  int status=0;
  fits_close_file(fptr, &status);
  // Don't throw if unwinding another exception:
  if (status && !std::uncaught_exception()) 
    throw_CFITSIO("closeFile() on " + getFilename());
}

template <class T>
void addColumn(FTable ft, string name, long repeat) {
  if (repeat<0 || repeat>1) ft.addColumn(vector<vector<T> >(), name, repeat);
  else ft.addColumn(vector<T>(), name, repeat);
}

template <class T>
void 
FITSTable::getFitsColumnData(FTable ft, string colName, int icol, 
			     long rowStart, long nRows) {
  if (nRows < 1) return;
  FITS::DataType dtype = FITS::FITSTypeOf<T>();
  int status=0;
  long repeat = ft.repeat(colName);
  if (repeat < 0) {
    // Branch for variable-length array
    // Get the vector that says how long each entry is
    long repeats[nRows];
    long offsets[nRows];
    // Note adding 1 to row numbers when FITS sees them, since CFITSIO is 1-indexed:
    fits_read_descripts(fptr, icol, rowStart+1, nRows, repeats, offsets, &status);
    long totalElements = 0;
    for (int i=0; i<nRows; i++)
      totalElements += repeats[i];
    vector<T> data(totalElements);
    T nullValue;
    int anyNulls;
    // Read all of the elements
    fits_read_col(fptr, dtype, icol, (LONGLONG) i+rowStart+1, (LONGLONG) 1, 
		  (LONGLONG) totalElements, &nullValue, &data[0], &anyNulls, &status);
    if (status) throw_CFITSIO("Reading table column " + colName);
    // Repackage into vectors of variable size
    vector<vector<T> > vv(nRows);
    typename vector<T>::const_iterator inptr = data.begin();
    for (int i=0; i<nRows; i++) {
      long length = repeats[i];
      vv[i] = vector<T>(inptr, inptr + length);
      inptr += length;
    }
    ft.writeCells(vv, colName, rowStart);
  } else if (repeat>1) {
    // Branch for fixed-length array
    long length=repeat * nRows;
    vector<T> data(length);
    T nullValue;
    int anyNulls;
    fits_read_col(fptr, dtype, icol, (LONGLONG) rowStart+1, (LONGLONG) 1, (LONGLONG) length,
		  &nullValue, &data[0], &anyNulls, &status);
    if (status) throw_CFITSIO("Reading table column " + colName);
    typename vector<T>::const_iterator inptr = data.begin();
    for (int i=0; i<nRows; i++) {
      ft.writeCell(vector<T>(inptr, inptr+repeat), colName, rowStart+i);
      inptr += repeat;
    }
  } else {
    // branch for scalar column
    T nullValue;
    int anyNulls;
    vector<T> data(nRows);
    fits_read_col(fptr, dtype, icol, (LONGLONG) rowStart+1, (LONGLONG) 1, (LONGLONG) nRows,
		  &nullValue, &data[0], &anyNulls, &status);
    if (status) throw_CFITSIO("Reading table column " + colName);
    ft.writeCells(data, colName, rowStart);
  }
}

FTable 
FITSTable::extract(long rowStart=0; long rowEnd=-1) const {
    img::FTable ft(nrows);
    vector<string> colNames;
    vector<int> colNums;
    moveTo();
    findColumns(ft, "*", colNames, colNums);
    createColumns(ft, colNames, colNums);
    readData(ft, colNames, colNums, rowStart, rowEnd);
    return ft;
}

FTable 
FITSTable::extract(const vector<string>& templates,
		   long rowStart=0; long rowEnd=-1) const {
    img::FTable ft(nrows);
    vector<string> colNames;
    vector<int> colNums;
    moveTo();
    for (int i=0; i<templates.size(); i++)
      findColumns(ft, templates[i], colNames, colNums);
    createColumns(ft, colNames, colNums);
    readData(ft, colNames, colNums, rowStart, rowEnd);
    return ft;
}

void
FITSTable::findColumns(FITSTable& ft,
		       string template,
		       vector<string>& names,
		       vector<int>& numbers) const {
    char colName[FITS::MAX_COLUMN_NAME_LENGTH];
    do {
      int colNum;
      fits_get_colname(fptr, CASEINSEN, template.c_str(), colName, &colNum, &status);
      if (status==0 || status==COL_NOT_UNIQUE) {
	// duplicate?
	for (int i=0; i<numbers.size(); i++)
	  if (colNum==numbers[i]) continue;
	// No; add to the lists
	colNames.push_back(colName);
	colNums.push_back(colNum);
      }
    } while (status==COL_NOT_UNIQUE);
    if (status==COL_NOT_FOUND) {
      status=0;
      fits_clear_errmsg();
    }
    if (status) throw_CFITSIO();
}

void
FITSTable::createColumns(FITSTable& ft, 
			 const vector<string>& names,
			 const vector<int>& numbers) const {
  int status=0;
  int datatype;
  long repeat, width;
  for (int i=0; i<colNames.size(); i++) {
    if (fits_get_eqcoltype(fptr, colNums[i], &datatype, &repeat, &width, &status))
      throw_CFITSIO("Getting info for column " + colName[i]);
    // width will give length of strings if this is string type
    bool isVariable = false;
    if (datatype < 0) {
      isVariable = true;
      datatype = -datatype;
    }
    
    // Add column to our table: big switch to decide which C type we will use
    // given what datatype says the FITS column data is
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
	// int is 4 bytes, that's what we want
	addColumn<int>(ft, colNames[i], repeat);
      else if (!FITS::CLongIsFITSLongLong)
	// If long is not 8 bytes, it's 4 and we'll use it
	addColumn<long>(ft, colNames[i], repeat);
      else
	throw FITSError("No intrinsic type holds 32-bit integers!");
      break;
      // Does not appear that unsigned int can arise, but put it here in case:
    case FITS::Tuint:
      addColumn<unsigned int>(ft, colNames[i], repeat);
      break;
    case FITS::Tulong:
      // ?? ambiguous what this would mean if it arose, go with intrinsic long
      addColumn<unsigned long>(ft, colNames[i], repeat);
      break;
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
      FormatAndThrow<FITSError>() << "Cannot convert FITS column data type " 
				  << dt << " to intrinsic type";
    } // end switch (dt)
  } // end loop over columns
  // ??? reserve ???
}

void 
FITSTable::readData(FITSTable& ft, 
		    const vector<string>& names,
		    const vector<int>& numbers,
		    long rowStart,
		    long rowEnd) const {

    // Now read in the data from FITS file into our columns, using the
    // recommended row increments for buffering
  int status = 0;
  long bufferRows;
  if (fits_get_rowsize(fitsptr, &bufferRows, &status))
    throw_CFITSIO();
  
  if (rowEnd < 0) rowEnd = nrows;
  for (long firstRow = rowStart; firstRow < rowEnd; firstRow+=bufferRows) {
    long rowCount = std::min(rowEnd, firstRow+bufferRows) - firstRow;
    for (int i=0; i<numbers.size(); i++) {
      string name = names[i];
      int num = numbers[i];

      switch (ft.elementType(name)) {
      case FITS::Tstring:
	GetFitsColumnData<string>(fitsptr, ft, name, num, rowStart, rowCount);
	break;
      case FITS::Tlogical:
	GetFitsColumnData<bool>(fitsptr, ft, name, num, rowStart, rowCount);
	break;
      case FITS::Tbyte:
	GetFitsColumnData<unsigned char>(fitsptr, ft, name, num, rowStart, rowCount);
	break;
      case FITS::Tsbyte:
	GetFitsColumnData<signed char>(fitsptr, ft, name, num, rowStart, rowCount);
	break;
      case FITS::Tshort:
	GetFitsColumnData<short>(fitsptr, ft, name, num, rowStart, rowCount);
	break;
      case FITS::Tushort:
	GetFitsColumnData<unsigned short>(fitsptr, ft, name, num, rowStart, rowCount);
	break;
      case FITS::Tint:
	GetFitsColumnData<int>(fitsptr, ft, name, num, rowStart, rowCount);
	break;
      case FITS::Tuint:
	GetFitsColumnData<unsigned int>(fitsptr, ft, name, num, rowStart, rowCount);
	break;
      case FITS::Tlong:
	GetFitsColumnData<long>(fitsptr, ft, name, num, rowStart, rowCount);
	break;
      case FITS::Tulong:
	GetFitsColumnData<unsigned long>(fitsptr, ft, name, num, rowStart, rowCount);
	break;
      case FITS::Tlonglong:
	GetFitsColumnData<LONGLONG>(fitsptr, ft, name, num, rowStart, rowCount);
	break;
      case FITS::Tfloat:
	GetFitsColumnData<float>(fitsptr, ft, name, num, rowStart, rowCount);
	break;
      case FITS::Tdouble:
	GetFitsColumnData<double>(fitsptr, ft, name, num, rowStart, rowCount);
	break;
      case FITS::Tcomplex:
	GetFitsColumnData<complex<float> >(fitsptr, ft, name, num, rowStart, rowCount);
	break;
      case FITS::Tdblcomplex:
	GetFitsColumnData<complex<double> >(fitsptr, ft, name, num, rowStart, rowCount);
	break;
      default:
	FormatAndThrow<FITSError>() << "Cannot convert data type " 
				    << ft.elementType(name) << " to FTable Column";
      }  // end switch on data type

    } // End column loop

  } // End row loop
