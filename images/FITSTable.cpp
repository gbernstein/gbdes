// Read/Writewrite FTables from FITS extensions

#include "FITSTable.h"

using namespace img;
using namespace FITS;

FITSTable::FITSTable(string filename, 
		     FITS::Flags f,
		     int hduNumber_): parent(filename, f), 
				      hduNumber(hduNumber_),
				      writeable( !(f==ReadOnly)) 
{
  if (hduNumber==0) 
    throw FITSError("Cannot open a FITS table at HDU number 0; file " + getFilename());
  int status=0;
  int hdutype;
  int nHDUs=0;
  // Count existing HDUs in this file
  // OK for empty file ???
  fits_get_num_hdus(parent, &nHDUs, &status);
  if (status) throw_CFITSIO();

  // Move to the requested HDU (request of -1 means to append new HDU)
  bool foundHDU = false;
  if (hduNumber>=0 && hduNumber < nHDUs) {
    // Try to find requested HDU, if one is requested
    // NOTE THAT CFITSIO HDU's are 1-indexed, I'm doing 0-indexed
    fits_movabs_hdu(parent, hduNumber+1, &hdutype, &status);
    if (status==BAD_HDU_NUM) {
      status = 0;
      fits_clear_errmsg();
    } else {
      foundHDU = true;
    }
  }
  if ( !foundHDU && (f & FITS::Create)) {
    // Try to make a new extension at the end of all of them or in requested spot

    // If we've asked for appending a new extension with hduNumber=-1,
    // then target an HDU number one past end
    if (hduNumber < 1) hduNumber = std::max(1,nHDUs);

    // Fill in 0-dimensional dummy to place our HDU in desired place.
    for (int emptyHDU=nHDUs; emptyHDU<hduNumber; emptyHDU++) {
      long naxes[2];
      fits_create_img(parent, BYTE_IMG, 0, naxes, &status);
    }

    // Now create new BINTABLE extension, empty initially
    {
      char* tDummy[1];
      tDummy[0] = 0;
      string hduName = "Change???";
      fits_create_tbl(parent, BINARY_TBL, (LONGLONG) 0, 0, 
		      tDummy, tDummy, tDummy,
		      const_cast<char*> (hduName.c_str()),
		      &status);
      // Make sure we can move to this
      fits_movabs_hdu(parent, hduNumber+1, &hdutype, &status);
      if (!status) foundHDU = true;
    }
  } else if (!foundHDU) {
    // Did not find HDU and not allowed to create it:
    FormatAndThrow<FITSError>() << "Did not find HDU " << hduNumber
				<< " in FITS file " << getFilename();
  } else if (f & FITS::OverwriteHDU) {
    // We did find our HDU number but are instructed to overwrite it,
    // so clean it out and put in a new empty bintable
    int dummy;
    fits_delete_hdu(parent, &dummy, &status);
    fits_flush_file(parent, &status);
    // ??? Note complications if we just deleted the primary header

    // Now insert new BINTABLE extension, empty initially, in desired spot
    // Need to move to the place *before* we want ours:
    // ???? Have some kind of bug (CFITSIO?) with deleting and replacing in the
    // middle of file; leaves something that has extra stuff at end of FITS file ????
    fits_movabs_hdu(parent, hduNumber, &hdutype, &status);
    {
      char* tDummy[1];
      tDummy[0] = 0;
      string hduName = "Change???";
      fits_insert_btbl(parent, (LONGLONG) 0, 0, 
		       tDummy, tDummy, tDummy,
		       const_cast<char*> (hduName.c_str()),
		       (long) 0,	// pcount: not reserving any space for heap.
		       &status);
      // Make sure we can move to this
      fits_movabs_hdu(parent, hduNumber+1, &hdutype, &status);
    }
  }
  // Report any errors getting here
  if (status) throw_CFITSIO("Error opening proper HDU for FITSTable file " +getFilename());
  // We should now have an opened HDU at desired number. See if it's a table as desired:
  if ( !(HDUType(hdutype) == HDUBinTable  ||
	 HDUType(hdutype) == HDUAsciiTable) )
    FormatAndThrow<FITSError>() << "HDU " << hduNumber
				<< " of file " << getFilename()
				<< " is not a FITS table extension";

  // Get the info that we want about this table:

  fits_get_num_cols(parent, &ncols, &status);
  fits_get_num_rows(parent, &nrows, &status);
  if (status) throw_CFITSIO("Opening FITSTable in file " + getFilename());
}

int
FITSTable::moveTo() const {
  int status=0;
  int hdutype;
  fits_movabs_hdu(parent, hduNumber+1, &hdutype, &status);
  return status;
}

FITSTable::~FITSTable() {
  // Flush this HDU if writeable, altered?
}

template <class T>
void addEmptyColumn(FTable ft, string name, long repeat, long stringLength=-1) {
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
      fits_get_colname(parent, CASEINSEN, 
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
      fits_get_eqcoltype(parent, colNumbers[i], &datatype, &repeat, &width, &status);
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
    /**cerr << "createColumn for <" << colNames[i] 
	     << "> column number " << colNumbers[i]
	     << " with dtype " << dt 
	     << " repeat " << repeat << " width " << width << endl; /***/
    switch (dt) {
    case FITS::Tbit:
      // I will convert bit to boolean on input
      addEmptyColumn<bool>(ft, colNames[i], repeat);
      break;
    case FITS::Tlogical:
      addEmptyColumn<bool>(ft, colNames[i], repeat);
      break;
    case FITS::Tbyte:
      addEmptyColumn<unsigned char>(ft, colNames[i], repeat);
      break;
    case FITS::Tsbyte:
      addEmptyColumn<signed char>(ft, colNames[i], repeat);
      break;
    case FITS::Tshort:
      addEmptyColumn<short>(ft, colNames[i], repeat);
      break;
    case FITS::Tushort:
      addEmptyColumn<unsigned short>(ft, colNames[i], repeat);
      break;
    case FITS::Tint32bit:
      if (FITS::CIntIsFITSLong)
	// int is 4 bytes, that's what we want
	addEmptyColumn<int>(ft, colNames[i], repeat);
      else if (!FITS::CLongIsFITSLongLong)
	// If long is not 8 bytes, it's 4 and we'll use it
	addEmptyColumn<long>(ft, colNames[i], repeat);
      else
	throw FITSError("No intrinsic type holds 32-bit integers!");
      break;
      // Does not appear that unsigned int can be in FITS Table, but put it here in case:
    case FITS::Tuint:
      addEmptyColumn<unsigned int>(ft, colNames[i], repeat);
      break;
    case FITS::Tulong:
      // ?? ambiguous what this would mean if it arose, go with intrinsic long
      addEmptyColumn<unsigned long>(ft, colNames[i], repeat);
      break;
    case FITS::Tlonglong:
      addEmptyColumn<LONGLONG>(ft, colNames[i], repeat);
      break;
    case FITS::Tfloat:
      addEmptyColumn<float>(ft, colNames[i], repeat);
      break;
    case FITS::Tdouble:
      addEmptyColumn<double>(ft, colNames[i], repeat);
      break;
    case FITS::Tcomplex:
      addEmptyColumn<complex<float> >(ft, colNames[i], repeat);
      break;
    case FITS::Tdblcomplex:
      addEmptyColumn<complex<double> >(ft, colNames[i], repeat);
      break;
    case FITS::Tstring:
      // Note that CFITSIO convention is that an array of fixed-length strings
      // can be specifed, in which case each string is length =width and repeat
      // gives total number of chars in whole array
      // If this was a variable-length char column, then I will call this a 
      // column with a single variable-length string in it:
      if (isVariable) addEmptyColumn<string>(ft, colNames[i], 1, -1);
      else if (width==repeat) {
	// Each cell is a single fixed-length string
	Assert(repeat>=1);
	addEmptyColumn<string>(ft, colNames[i], 1, repeat);
      } else {
	// Each cell is a fixed-size array of fixed-length strings
	int nStrings = repeat / width;
	Assert(nStrings>1);
	addEmptyColumn<string>(ft, colNames[i], nStrings, width);
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
    fits_get_rowsize(parent, &bufferRows, &status);
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
				    << ft.elementType(name) << " to FTable for column "
				    << name;
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
	fits_read_descript(parent, icol, rowStart+1+i, &nElements, &offset, &status);
	data.resize(nElements);
	fits_read_col(parent, dtype, icol, (LONGLONG) rowStart+1+i, (LONGLONG) 1, 
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
      fits_read_col(parent, dtype, icol, (LONGLONG) rowStart+1, (LONGLONG) 1, (LONGLONG) nElements,
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
      fits_read_col(parent, dtype, icol, (LONGLONG) rowStart+1, (LONGLONG) 1, (LONGLONG) nRows,
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
	fits_read_descript(parent, icol, rowStart+1+i, &nElements, &offset, &status);
	data.resize(nElements);
	fits_read_col(parent, dtype, icol, (LONGLONG) rowStart+1+i, (LONGLONG) 1, 
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
      fits_read_col(parent, dtype, icol, (LONGLONG) rowStart+1, (LONGLONG) 1, (LONGLONG) nElements,
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
      fits_read_col(parent, dtype, icol, (LONGLONG) rowStart+1, (LONGLONG) 1, (LONGLONG) nRows,
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
	  fits_read_descript(parent, icol, rowStart+1+i, &nElements, &offset, &status);
	  data = new char[nElements+1];
	  // Read as a c-string; wants char** for destination
	  fits_read_col(parent, Tbyte, icol, (LONGLONG) rowStart+1+i, (LONGLONG) 1, 
			(LONGLONG) nElements, nullPtr, &data, &anyNulls, &status);
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
      fits_read_col(parent, dtype, icol, (LONGLONG) rowStart+1, (LONGLONG) 1, (LONGLONG) nStrings,
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


////////////////////////////////////////////////////////////////////////////////
////////// Writing routines
////////////////////////////////////////////////////////////////////////////////

// Write contents of FTable into the FITSTable
void
FITSTable::replaceWith(FTable ft) {
  clear();
  long rowStart = 0;
  long rowEnd = ft.nrows();
  vector<string> colNames;
  vector<int> colNums;
  addFitsColumns(ft, colNames, colNums);
  // Now going to write out the contents to FITS, using recommended buffering row intervals.
  int status = 0;
  long bufferRows;
#pragma omp critical (fitsio)
  {
    status = moveTo();
    fits_get_rowsize(parent, &bufferRows, &status);
  }
  if (status) throw_CFITSIO();
  
  for (long firstRow = rowStart; firstRow < rowEnd; firstRow+=bufferRows) {
    long rowCount = std::min(rowEnd, firstRow+bufferRows) - firstRow;
    // Loop over columns, writing interval of each to output
    for (int iCol = 0; iCol < colNames.size(); iCol++) {
      FITS::DataType dt = ft.elementType(colNames[iCol]);
      // Dispatch the data according to the datatype:
      switch (dt) {
      case Tlogical:
	writeFitsColumn<bool>(ft, colNames[iCol], colNums[iCol], 
			      firstRow, firstRow+rowCount);
	break;
      case Tstring:
	writeFitsColumn<string>(ft, colNames[iCol], colNums[iCol], 
				firstRow, firstRow+rowCount);
	break;
      case Tbyte:
	writeFitsColumn<unsigned char>(ft, colNames[iCol], colNums[iCol], 
				       firstRow, firstRow+rowCount);
	break;
      case Tsbyte:
	writeFitsColumn<signed char>(ft, colNames[iCol], colNums[iCol], 
				     firstRow, firstRow+rowCount);
	break;
      case Tshort:
	writeFitsColumn<short>(ft, colNames[iCol], colNums[iCol], 
			       firstRow, firstRow+rowCount);
	break;
      case Tushort:
	writeFitsColumn<unsigned short>(ft, colNames[iCol], colNums[iCol], 
					firstRow, firstRow+rowCount);
	break;
      case Tint:
	writeFitsColumn<int>(ft, colNames[iCol], colNums[iCol], 
			     firstRow, firstRow+rowCount);
	break;
      case Tuint:
	writeFitsColumn<unsigned int>(ft, colNames[iCol], colNums[iCol], 
				      firstRow, firstRow+rowCount);
	break;
      case Tlong:
	writeFitsColumn<long>(ft, colNames[iCol], colNums[iCol], 
			     firstRow, firstRow+rowCount);
	break;
      case Tulong:
	writeFitsColumn<unsigned long>(ft, colNames[iCol], colNums[iCol], 
				      firstRow, firstRow+rowCount);
	break;
      case Tlonglong:
	writeFitsColumn<LONGLONG>(ft, colNames[iCol], colNums[iCol], firstRow, firstRow+rowCount);
	break;
      case Tfloat:
	writeFitsColumn<float>(ft, colNames[iCol], colNums[iCol], firstRow, firstRow+rowCount);
	break;
      case Tdouble:
	writeFitsColumn<double>(ft, colNames[iCol], colNums[iCol], firstRow, firstRow+rowCount);
	break;
      case Tcomplex:
	writeFitsColumn<complex<float> >(ft, colNames[iCol], colNums[iCol], 
					 firstRow, firstRow+rowCount);
	break;
      case Tdblcomplex:
	writeFitsColumn<complex<double> >(ft, colNames[iCol], colNums[iCol], 
					  firstRow, firstRow+rowCount);
	break;
      default:
	FormatAndThrow<FITSError>() << "Unknown DataType " << dt
				    << " for writing data to FITS column";
      } // end DataType switch
    } // end column loop
  } // end row chunk loop
}

void
FITSTable::addFitsColumns(FTable ft, vector<string>& colNames, vector<int>& colNums) {
  if (!isWriteable())
    throw FITSError("Attempt to add columns to read-only HDU in FITS file " + getFilename());
  int status=0;
  int colNum;
  // Add these columns to end of any there:
#pragma omp critical (fitsio)
  {
    status = moveTo();
    fits_get_num_cols(parent, &colNum, &status);
  }
  if (status) throw_CFITSIO();
  colNum++;	// columns are 1-indexed, and this arg gives desired posn of insertion
  colNums.clear();

  // Get all column names
  colNames = ft.columnNames();
  for (int i=0; i<colNames.size(); i++, colNum++) {
    // loop over ft's columns:
    // First make the TFORM string
    long repeat = ft.repeat(colNames[i]);
    ostringstream oss;
    if (ft.elementType(colNames[i]) == FITS::Tstring) {
      long length = ft.stringLength(colNames[i]);
      if (length < 0) {
	// single variable-length string per cell:
	if (repeat==1) oss << "1PA";
	// otherwise: we will punt right now
	throw FITSError("Cannot format variable-length string array for bintable "
			"at column " + colNames[i]);
      } else {
	// Fixed-length strings
	if (repeat < 0) {
	  throw FITSError("Cannot format variable-length arrays of strings "
			  " for bintable at column " + colNames[i]);
	} else if (repeat==1) {
	  // Single string
	  oss << length << "A";
	} else {
	  // fixed-length array.  Use CFITSIO convention
	  oss << length*repeat << "A" << length;
	}
      }
    } else {
      // For all non-string types:
      if (repeat<0) oss << "1P";
      else oss << repeat;
      oss << ft.columnCode(colNames[i]);
    }

#pragma omp critical (fitsio)
    {
      status = moveTo();
      fits_insert_col(parent, colNum, const_cast<char*> (colNames[i].c_str()), 
		      const_cast<char*> (oss.str().c_str()),
		      &status);
    }
    colNums.push_back(colNum);
    if (status) throw_CFITSIO("creating column " + colNames[i]);
  }
}

void
FITSTable::clear() {
  // Remove all data from a FITSTable (no flush)
  // ??? clear header too ???
  if (!isWriteable())
    throw FITSError("Attempt to clear read-only HDU in FITS file " + getFilename());
  int status;
#pragma omp critical (fitsio)
  {
    status = moveTo();
    int ncols;
    // Delete all columns, from back to front:
    fits_get_num_cols(parent, &ncols, &status);
    for (int i=ncols; i>0; i--)
      fits_delete_col(parent, i, &status);
    // compress (=empty) the heap:
    fits_compress_heap(parent, &status);
  }
  if (status) throw_CFITSIO("In FITSTable::clear(): ");
}

template <class DT> 
void
FITSTable::writeFitsColumn(FTable ft, string colName, int colNum,
			   long rowStart, long rowEnd) {
  if (!isWriteable())
    throw FITSError("Attempt to write to read-only HDU in FITS file " + getFilename());
  int status=0;
  long repeat = ft.repeat(colName);
  vector<DT> data;
  if (repeat==1) {
    // Scalar column, can write whole group at once
    ft.readCells(data, colName, rowStart, rowEnd);
#pragma omp critical (fitsio)
    {
      status = moveTo();
      fits_write_col(parent, FITS::FITSTypeOf<DT>(), colNum, (LONGLONG) rowStart+1,
		     (LONGLONG) 1, (LONGLONG) (rowEnd-rowStart), &data[0], &status);
    }
    if (status) throw_CFITSIO("Writing to FITS column " + colName);
  } else {
    // Array cells, write 1 row at a time
    for (int iRow = rowStart; iRow < rowEnd; iRow++) {
      ft.readCell(data, colName, iRow);
      if (repeat >=0) Assert(data.size()==repeat);
#pragma omp critical (fitsio)
      {
	status = moveTo();
	fits_write_col(parent, FITS::FITSTypeOf<DT>(), colNum, (LONGLONG) iRow+1,
		       (LONGLONG) 1, (LONGLONG) data.size(), &data[0], &status);
      }
      if (status) throw_CFITSIO("Writing to FITS column " + colName);
    } // end row loop
  } // End if/else for array/scalar cells
}

// Specialize for logical, need to pass data in/out as char
template <> 
void
FITSTable::writeFitsColumn<bool>(FTable ft, string colName, int colNum,
				 long rowStart, long rowEnd) {
  if (!isWriteable())
    throw FITSError("Attempt to write to read-only HDU in FITS file " + getFilename());
  int status=0;
  long repeat = ft.repeat(colName);
  vector<bool> data;
  if (repeat==1) {
    // Scalar column, can write whole group at once
    ft.readCells(data, colName, rowStart, rowEnd);
    char cdata[data.size()];
    for (int i=0; i<data.size(); i++)
      cdata[i] = data[i] ? 1 : 0;
#pragma omp critical (fitsio)
    {
      status = moveTo();
      fits_write_col(parent, FITS::Tlogical, colNum, (LONGLONG) rowStart+1,
		     (LONGLONG) 1, (LONGLONG) (rowEnd-rowStart), cdata, &status);
    }
    if (status) throw_CFITSIO("Writing to FITS column " + colName);
  } else {
    // Array cells, write 1 row at a time
    for (int iRow = rowStart; iRow < rowEnd; iRow++) {
      ft.readCell(data, colName, iRow);
      char cdata[data.size()];
      for (int i=0; i<data.size(); i++)
	cdata[i] = data[i] ? 1 : 0;
      if (repeat >=0) Assert(data.size()==repeat);
#pragma omp critical (fitsio)
      {
	status = moveTo();
	fits_write_col(parent, FITS::Tlogical, colNum, (LONGLONG) iRow+1,
		       (LONGLONG) 1, (LONGLONG) data.size(), cdata, &status);
      }
      if (status) throw_CFITSIO("Writing to FITS column " + colName);
    } // end row loop
  } // End if/else for array/scalar cells
}

// Specialize for strings
template<>
void
FITSTable::writeFitsColumn<string>(FTable ft, string colName, int colNum,
					long rowStart, long rowEnd) {
  if (!isWriteable())
    throw FITSError("Attempt to write to read-only HDU in FITS file " + getFilename());
  int status=0;
  long repeat = ft.repeat(colName);
  long length = ft.stringLength(colName);
  if (repeat==1) {
    // One string per cell.  Read individually, whether or not fixed-length: FITSIO will deal
    string data;
    for (int iRow = rowStart; iRow < rowEnd; iRow++) {
      ft.readCell(data, colName, iRow);
      char* cdata[1];
      cdata[0] = const_cast<char*> (data.c_str());
#pragma omp critical (fitsio)
      {
	status = moveTo();
	// fits_write_col wants char** input for source of data
	fits_write_col(parent, Tstring, colNum, (LONGLONG) iRow+1,
		       (LONGLONG) 1, (LONGLONG) 1,
		       cdata, &status);
      }
      if (status) throw_CFITSIO("Writing to FITS column " + colName);
    }
  } else {
      if (length < 0 || repeat < 0)
	throw FITSError("Cannot write column <" + colName + "> to FITS because"
			" it has variable-length string arrays.");
      // Array cells, write 1 row at a time
      vector<string> data;
      char* cdata[repeat];
      for (int iRow = rowStart; iRow < rowEnd; iRow++) {
	ft.readCell(data, colName, iRow);
	data.resize(repeat);
	for (int j=0; j<repeat; j++)
	  cdata[j] = const_cast<char*> (data[j].c_str());
#pragma omp critical (fitsio)
      {
	status = moveTo();
	fits_write_col(parent, Tstring, colNum, (LONGLONG) iRow+1,
		       (LONGLONG) 1, (LONGLONG) repeat, cdata, &status);
      }
      if (status) throw_CFITSIO("Writing to FITS column " + colName);
    } // end row loop
  } // End if/else for array/scalar cells
}



