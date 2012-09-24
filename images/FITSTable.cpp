// Read/Writewrite FTables from FITS extensions

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
  if ( !(HDUType(hdutype) == HDUBinTable  ||
       HDUType(hdutype) == HDUAsciiTable) )
    throw FITSError("Not a FITS table extension");

  fits_get_num_cols(fptr, &ncols, &status);
  fits_get_num_rows(fptr, &nrows, &status);
  if (status) throw_CFITSIO("Opening FITSTable in file " + fname);
}

int
FITSTable::moveTo() const {
  int status=0;
  int hdutype;
  fits_movabs_hdu(fptr, hduNumber, &hdutype, &status);
  return status;
}

FITSTable::~FITSTable() {
  if (!fptr) return;
  int status=0;
#pragma openmp critical (fitsio)
  {
    fits_close_file(fptr, &status);
  }
  // Don't throw if unwinding another exception:
  if (status && !std::uncaught_exception()) 
    throw_CFITSIO("closeFile() on " + getFilename());
}

template <class T>
void addColumn(FTable ft, string name, long repeat, long stringLength=-1) {
  if (repeat<0 || repeat>1) ft.addColumn(vector<vector<T> >(), name, repeat, stringLength);
  else ft.addColumn(vector<T>(), name, repeat, stringLength);
}


FTable 
FITSTable::extract(long rowStart, long rowEnd) const {
  if (rowEnd < 0 ) rowEnd = nrows;
  if (rowEnd < rowStart) rowEnd = rowStart;  // Makes an empty table
  int outRows = rowEnd - rowStart;
  img::FTable ft(outRows);
  vector<string> colNames;
  vector<int> colNums;
  const string matchAll("*");
  findColumns(ft, matchAll, colNames, colNums);
  createColumns(ft, colNames, colNums);
  readData(ft, colNames, colNums, rowStart, rowEnd);
  return ft;
}

FTable 
FITSTable::extract(const vector<string>& templates,
		   long rowStart, long rowEnd) const {
  if (rowEnd < 0 ) rowEnd = nrows;
  if (rowEnd < rowStart) rowEnd = rowStart;  // Makes an empty table
  int outRows = rowEnd - rowStart;
  img::FTable ft(outRows);
  vector<string> colNames;
  vector<int> colNums;
  for (int i=0; i<templates.size(); i++)
    findColumns(ft, templates[i], colNames, colNums);
  createColumns(ft, colNames, colNums);
  readData(ft, colNames, colNums, rowStart, rowEnd);
  return ft;
}

void
FITSTable::findColumns(FTable& ft,
		       string matchMe,
		       vector<string>& colNames,
		       vector<int>& colNumbers) const {
  char colName[FITS::MAX_COLUMN_NAME_LENGTH];

  int status = 0;
#pragma omp critical (fitsio) 
  {  
    status = moveTo();
    do {
      int colNum;
      // Note const_cast to match non-const template for CFITSIO
      fits_get_colname(fptr, CASEINSEN, 
		       const_cast<char*> (matchMe.c_str()), colName, &colNum, &status);
      if (status==0 || status==COL_NOT_UNIQUE) {
	// duplicate?
	for (int i=0; i<colNumbers.size(); i++)
	  if (colNum==colNumbers[i]) continue;
	// No; add to the lists
	colNames.push_back(colName);
	colNumbers.push_back(colNum);
      }
    } while (status==COL_NOT_UNIQUE);
    if (status==COL_NOT_FOUND) {
      status=0;
      fits_clear_errmsg();
    }
  } // end omp critical region
  if (status) throw_CFITSIO();
}

void
FITSTable::createColumns(FTable& ft, 
			 const vector<string>& colNames,
			 const vector<int>& colNumbers) const {
  int status=0;
  int datatype;
  long repeat, width;
  for (int i=0; i<colNames.size(); i++) {
#pragma omp critical (fitsio) 
    {  
      status = moveTo();
      fits_get_eqcoltype(fptr, colNumbers[i], &datatype, &repeat, &width, &status);
    }
    if (status) throw_CFITSIO("Getting info for column " + colNames[i]);

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
    case FITS::Tbit:
      // I will convert bit to boolean on input
      addColumn<bool>(ft, colNames[i], repeat);
      break;
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
      // Does not appear that unsigned int can be in FITS Table, but put it here in case:
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
      // Note that CFITSIO convention is that an array of fixed-length strings
      // can be specifed, in which case each string is length =width and repeat
      // gives total number of chars in whole array
      // If this was a variable-length char column, then I will call this a 
      // column with a single variable-length string in it:
      if (isVariable) addColumn<string>(ft, colNames[i], 1, -1);
      else if (width==repeat) {
	// Each cell is a single fixed-length string
	Assert(repeat>=1);
	addColumn<string>(ft, colNames[i], 1, repeat);
      } else {
	// Each cell is a fixed-size array of fixed-length strings
	int nStrings = repeat / width;
	Assert(nStrings>1);
	addColumn<vector<string> >(ft, colNames[i], nStrings, width);
      }
      break;
    default:
      FormatAndThrow<FITSError>() << "Cannot convert FITS column data type " 
				  << dt << " to intrinsic type";
    } // end switch (dt)
  } // end loop over columns
}

void 
FITSTable::readData(FTable& ft, 
		    const vector<string>& names,
		    const vector<int>& numbers,
		    long rowStart,
		    long rowEnd) const {

  // Now read in the data from FITS file into our columns, using the
  // recommended row increments for buffering
  int status = 0;
  long bufferRows;
#pragma omp critical (fitsio)
  {
    status = moveTo();
    fits_get_rowsize(fptr, &bufferRows, &status);
  }
  if (status) throw_CFITSIO();
  
  for (long firstRow = rowStart; firstRow < rowEnd; firstRow+=bufferRows) {
    long rowCount = std::min(rowEnd, firstRow+bufferRows) - firstRow;
    for (int i=0; i<numbers.size(); i++) {
      string name = names[i];
      int num = numbers[i];

      switch (ft.elementType(name)) {
      case FITS::Tstring:
	getFitsColumnData<string>(ft, name, num, rowStart, rowCount);
	break;
      case FITS::Tlogical:
	getFitsColumnData<bool>(ft, name, num, rowStart, rowCount);
	break;
      case FITS::Tbyte:
	getFitsColumnData<unsigned char>(ft, name, num, rowStart, rowCount);
	break;
      case FITS::Tsbyte:
	getFitsColumnData<signed char>(ft, name, num, rowStart, rowCount);
	break;
      case FITS::Tshort:
	getFitsColumnData<short>(ft, name, num, rowStart, rowCount);
	break;
      case FITS::Tushort:
	getFitsColumnData<unsigned short>(ft, name, num, rowStart, rowCount);
	break;
      case FITS::Tint:
	getFitsColumnData<int>(ft, name, num, rowStart, rowCount);
	break;
      case FITS::Tuint:
	getFitsColumnData<unsigned int>(ft, name, num, rowStart, rowCount);
	break;
      case FITS::Tlong:
	getFitsColumnData<long>(ft, name, num, rowStart, rowCount);
	break;
      case FITS::Tulong:
	getFitsColumnData<unsigned long>(ft, name, num, rowStart, rowCount);
	break;
      case FITS::Tlonglong:
	getFitsColumnData<LONGLONG>(ft, name, num, rowStart, rowCount);
	break;
      case FITS::Tfloat:
	getFitsColumnData<float>(ft, name, num, rowStart, rowCount);
	break;
      case FITS::Tdouble:
	getFitsColumnData<double>(ft, name, num, rowStart, rowCount);
	break;
      case FITS::Tcomplex:
	getFitsColumnData<complex<float> >(ft, name, num, rowStart, rowCount);
	break;
      case FITS::Tdblcomplex:
	getFitsColumnData<complex<double> >(ft, name, num, rowStart, rowCount);
	break;
      default:
	FormatAndThrow<FITSError>() << "Cannot convert data type " 
				    << ft.elementType(name) << " to FTable Column";
      }  // end switch on data type

    } // End column loop

  } // End row loop

} // end readData()

// Now a template for routine that gets our column data from file via CFITSIO routines
// Note that this is being called with template being the type of the data in the
// C++ FTable.
template <class T>
void 
FITSTable::getFitsColumnData(FTable ft, string colName, int icol, 
			     long rowStart, long nRows) const {
  if (nRows < 1) return;
  FITS::DataType dtype = FITS::FITSTypeOf<T>();
  int status=0;
  T* nullPtr = 0;
  int anyNulls;
  long repeat = ft.repeat(colName);
  if (repeat < 0) {
    // Branch for variable-length array: read row by row
    long nElements;
    long offset;
    vector<T> data;
    for (int i=0; i<nRows; i++) {
#pragma omp critical (fitsio)
      {
	status = moveTo();
	// Note adding 1 to row numbers when FITS sees them, since CFITSIO is 1-indexed:
	fits_read_descript(fptr, icol, rowStart+1+i, &nElements, &offset, &status);
	data.resize(nElements);
	fits_read_col(fptr, dtype, icol, (LONGLONG) rowStart+1+i, (LONGLONG) 1, 
		      (LONGLONG) nElements, nullPtr, &data[0], &anyNulls, &status);
      }
      if (status) throw_CFITSIO("Reading table column " + colName);
      ft.writeCell(data, colName, rowStart+i);
    } // end row loop
  } else if (repeat>1) {
    // Branch for fixed-length array - can read all rows at one time
    long nElements = repeat * nRows;
    vector<T> data(nElements);
#pragma omp critical (fitsio)
    {
      status = moveTo();
      fits_read_col(fptr, dtype, icol, (LONGLONG) rowStart+1, (LONGLONG) 1, (LONGLONG) nElements,
		    nullPtr, &data[0], &anyNulls, &status);
    }
    if (status) throw_CFITSIO("Reading table column " + colName);
    typename vector<T>::const_iterator inptr = data.begin();
    for (int i=0; i<nRows; i++) {
      ft.writeCell(vector<T>(inptr, inptr+repeat), colName, rowStart+i);
      inptr += repeat;
    }
  } else {
    // branch for scalar column
    vector<T> data(nRows);
#pragma omp critical (fitsio)
    {
      status = moveTo();
      fits_read_col(fptr, dtype, icol, (LONGLONG) rowStart+1, (LONGLONG) 1, (LONGLONG) nRows,
		    nullPtr, &data[0], &anyNulls, &status);
    }
    if (status) throw_CFITSIO("Reading table column " + colName);
    ft.writeCells(data, colName, rowStart);
  }
}


// Specialize for logical: CFITSIO is expecting char arrays to hold data, so 
// just do this explicitly instead of assuming that bool = char on this system.
template <>
void 
FITSTable::getFitsColumnData<bool>(FTable ft, string colName, int icol, 
				   long rowStart, long nRows) const {
  if (nRows < 1) return;
  FITS::DataType dtype = FITS::Tlogical;
  int status=0;
  char* nullPtr = 0;
  int anyNulls;
  long repeat = ft.repeat(colName);
  if (repeat < 0) {
    // Branch for variable-length array: read row by row
    long nElements;
    long offset;
    vector<char> data;
    for (int i=0; i<nRows; i++) {
#pragma omp critical (fitsio)
      {
	status = moveTo();
	// Note adding 1 to row numbers when FITS sees them, since CFITSIO is 1-indexed:
	fits_read_descript(fptr, icol, rowStart+1+i, &nElements, &offset, &status);
	data.resize(nElements);
	fits_read_col(fptr, dtype, icol, (LONGLONG) rowStart+1+i, (LONGLONG) 1, 
		      (LONGLONG) nElements, nullPtr, &data[0], &anyNulls, &status);
      }
      if (status) throw_CFITSIO("Reading table column " + colName);
      // Rewrite as bool
      ft.writeCell(vector<bool>(data.begin(), data.end()), colName, rowStart+i);
    } // end row loop
  } else if (repeat>1) {
    // Branch for fixed-length array
    long nElements=repeat * nRows;
    vector<char> data(nElements);
#pragma omp critical (fitsio)
    {
      status = moveTo();
      fits_read_col(fptr, dtype, icol, (LONGLONG) rowStart+1, (LONGLONG) 1, (LONGLONG) nElements,
		    nullPtr, &data[0], &anyNulls, &status);
    }
    if (status) throw_CFITSIO("Reading table column " + colName);
    typename vector<char>::const_iterator inptr = data.begin();
    for (int i=0; i<nRows; i++) {
      ft.writeCell(vector<bool>(inptr, inptr+repeat), colName, rowStart+i);
      inptr += repeat;
    }
  } else {
    // branch for scalar column
    vector<char> data(nRows);
#pragma omp critical (fitsio)
    {
      status = moveTo();
      fits_read_col(fptr, dtype, icol, (LONGLONG) rowStart+1, (LONGLONG) 1, (LONGLONG) nRows,
		    nullPtr, &data[0], &anyNulls, &status);
    }
    if (status) throw_CFITSIO("Reading table column " + colName);
    // Convert to a bool vector
    vector<bool> booldata(data.begin(), data.end());
    ft.writeCells(booldata, colName, rowStart);
  }
}

// Specialize for string
// Complications because bintables store strings as char arrays.  So for FITSIO, "repeat" is
// of chars, not number of strings, but FTable will instead use "stringLength()".  
// There is a hack (use TDIM) to store an array of fixed
// number of strings of fixed length (can fill with nulls to use smaller strings.
// Not sure whether variable-length char arrays are ever used.
template <>
void 
FITSTable::getFitsColumnData<string>(FTable ft, string colName, int icol, 
				     long rowStart, long nRows) const {
  if (nRows < 1) return;
  FITS::DataType dtype = Tstring;
  int status=0;
  char* nullPtr = 0;
  int anyNulls;
  long repeat = ft.repeat(colName);
  // First consider variable-length strings:
  long length = ft.stringLength(colName);
  if (length < 0) {
    if (repeat == 1) {
      // Read a single variable-length string for each row
      long nElements;
      long offset;
      char *data;
      for (int i=0; i<nRows; i++) {
#pragma omp critical (fitsio)
	{
	  status = moveTo();
	  // Note adding 1 to row numbers when FITS sees them, since CFITSIO is 1-indexed:
	  fits_read_descript(fptr, icol, rowStart+1+i, &nElements, &offset, &status);
	  data = new char[nElements+1];
	  // Read as a byte (char) array ???
	  fits_read_col(fptr, Tbyte, icol, (LONGLONG) rowStart+1+i, (LONGLONG) 1, 
			(LONGLONG) nElements, nullPtr, data, &anyNulls, &status);
	}
	if (status) throw_CFITSIO("Reading table column " + colName);
	// Rewrite as string: make sure it's null-terminated
	data[nElements] = 0;
	string s(data);
	ft.writeCell(s, colName, rowStart+i);
	delete[] data;
      } // end row loop
    } else {
    // Don't know how FITS would store arrays of variable-length strings!
    FormatAndThrow<FITSError>() << "FITS table retrieval of arrays of variable-length strings"
				<< " not implemented.  Column " << colName;
    }
  } else {
    // Branch for fixed-length string per cell
    long nStrings= repeat * nRows;
    vector<char*> data(nStrings);
    for (int i=0; i<nStrings; i++)
      data[i] = new char[length+1];
#pragma omp critical (fitsio)
    {
      status = moveTo();
      // ??? nStrings or # characters here ???
      fits_read_col(fptr, dtype, icol, (LONGLONG) rowStart+1, (LONGLONG) 1, (LONGLONG) nStrings,
		    nullPtr, &data[0], &anyNulls, &status);
    }
    if (status) throw_CFITSIO("Reading table column " + colName);

    if (repeat==1) {
      // One string per cell; make an array of them all
      vector<string> vs(nRows);
      for (int i=0; i<nRows; i++) {
	// Make sure it's null terminated:
	data[i][length]=0;
	vs[i] = data[i];
	delete[] data[i];
      }
      ft.writeCells(vs, colName, rowStart);
    } else {
      // Multiple strings per cell: write cells one at a time
      vector<string> vs(repeat);
      int iString = 0;
      for (int j=0; j<nRows; j++) {
	for (int i=0; i<repeat; i++, iString++) {
	  data[iString][length] = 0; // insure null termination
	  vs[i] = data[iString];
	  delete[] data[iString];
	}
	ft.writeCell(vs, colName, rowStart+j);
      }
    }
  } // end if/then/else for array configuration
}
