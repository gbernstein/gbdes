// Read/Writewrite FTables from FITS extensions

#include "FitsTable.h"
#include "Expressions.h"

using namespace img;
using namespace FITS;

FitsTable::FitsTable(string filename, 
		     FITS::Flags f,
		     int hduNumber_): Hdu(filename, HDUBinTable, hduNumber_, f),
				      dptr(0), dcount(0)  {}

FitsTable::FitsTable(string filename, 
		     FITS::Flags f,
		     string hduName_): Hdu(filename, HDUBinTable, hduName_, f),
				      dptr(0), dcount(0)  {}

FitsTable::~FitsTable() {
  if (dptr) {
    // Make sure nothing else is using the mirrored data still:
    if (*dcount != 1)
      throwFitsOrDump("Closing FitsTable with its data still linked, file " + getFilename());
    if (isWriteable() && dptr->isChanged()) flushData();
    delete dptr;
    delete dcount;
  }
  // Destructor for Hdu will flush header if still needed, and close files.
}

std::vector<std::string> 
FitsTable::listColumns() const {
  if (dptr)
    return dptr->listColumns();
  else 
    return listFitsColumns();
}

img::FTable
FitsTable::use() {
  if (!dptr) {
    dptr = loadData();
    dcount = new int(1);
  }
  Hdu::header();	// Force header loading
  // Build the table from our data and header:
  img::FTable ft(hptr, hcount, dptr, dcount);
  if (!isWriteable()) ft.setLock();
  return ft;
}

// Here is a dummy tokenizer that just finds names of all columns we
// will need data from.  Don't try to parse after tokenizing with this!
class FindColumnsToken: public expressions::Token {
public:
  FindColumnsToken(const Header* hptr_): hptr(hptr_) {}
  const set<string>& columnNames() const {return columns;}
  virtual expressions::Token* createFromString(const std::string& input,
					      size_t& begin, size_t& end,
					      bool lastTokenWasOperator) const {
    Assert(end>begin);
    if (input[begin]=='@') {
      // Read a header keyword as a constant
      string keyword = input.substr(begin+1, end-(begin+1));
      bool vb;
      string vs;
      int vi;
      double vd;
      if (hptr->getValue(keyword, vb)
	  || hptr->getValue(keyword, vs)
	  || hptr->getValue(keyword, vi)
	  || hptr->getValue(keyword, vd) )
	return new FindColumnsToken(hptr);
      else 
	throw FTableError("FitsTable expression evaluation cannot find Header keyword "
			  + keyword);
    } else {
      string colName = input.substr(begin,end-begin);
      columns.insert(colName);
      return new FindColumnsToken(hptr);
    }
  }
private:
  const Header* hptr;
  mutable std::set<std::string> columns;
};


img::FTable 
FitsTable::extract(const string& expression,
		   const std::vector<std::string>& colMatches) const { 
  img::TableData* tdata = new TableData();
  if (dptr) {
    vector<bool> vb;
    dptr->evaluate(vb, expression, header());
    dptr->filterRows(tdata, vb);
  } else {
    FindColumnsToken tokenReader(header());

    // tokenize just to find out what columns are needed
    std::list<expressions::Token*> 
      tokenList = expressions::tokenize(expression, tokenReader);

    // Get rid of these tokens, we'll re-parse below
    for (std::list<expressions::Token*>::iterator i = tokenList.begin();
	 i != tokenList.end();
	 ++i)
      delete *i; // ??? Need to do better...

    // Read from file all desired columns plus those needed for filtering
    const set<string>& forFiltering = tokenReader.columnNames();
    vector<string> augmented = colMatches;
    for ( set<string>::const_iterator i = forFiltering.begin();
	  i != forFiltering.end();
	  ++i)
      augmented.push_back(*i);

    TableData* tdata2 = loadData(0, -1, augmented);
    vector<bool> vb;
    tdata2->evaluate(vb, expression, header());
    tdata2->filterRows(tdata, vb);
    delete tdata2;
  }

  // Retain only the originally request columns
  tdata->filter(tdata, 0, -1, colMatches);
  img::Header* hh = header()->duplicate();
  return FTable(hh, new int(0), tdata, new int(0)); 
}

img::FTable 
FitsTable::extract(long rowStart, long rowEnd,
		   const std::vector<std::string>& colMatches) const { 
  img::TableData* tdata;
  if (dptr) {
    tdata = dptr->extract(rowStart, rowEnd, colMatches);
  } else {
    tdata = loadData(rowStart, rowEnd, colMatches);
  }
  img::Header* hh = header()->duplicate();
  return FTable(hh, new int(0), tdata, new int(0)); 
}

void
FitsTable::copy(FTable ft) {
  checkWriteable("FitsTable::copy()");
  if (hptr) {
    hptr->copyFrom(*ft.header());
  } else {
    hptr = ft.header()->duplicate();
    hcount = new int(1);
  }
  img::Header* hh;
  int* hc;
  img::TableData* dd;
  int* dc;
  ft.showGuts(hh, hc, dd, dc);
  if (dptr) {
    dptr->copyFrom(*dd);
  } else {
    dptr = dd->duplicate();
    dcount = new int(1);
  }
  // Mark both as altered
  hptr->touch();
  dptr->touch();
}
 
// Attach this FITS file to data & header from external FTable:
void
FitsTable::adopt(FTable ft) { 
  checkWriteable("FitsTable::adopt()");
  img::Header* hh;
  int* hc;
  img::TableData* dd;
  int* dc;
  ft.showGuts(hh, hc, dd, dc);
  // Attach to header using Hdu base class:
  adoptHeader(hh, hc);
  // Get rid of any old data:
  if (dcount) {
    (*dcount)--;
    if (*dcount<=0) {
      delete dptr;
      delete dcount;
    }
  }
  dptr = dd;
  dcount = dc;
  (*dcount)++;
  dptr->touch();
}

void
FitsTable::flushData() {
  // Flush both data to FITS extension:
  if (!isWriteable()) return;
  if (dptr && dptr->isChanged())
    writeFitsTableData(dptr);
  dptr->clearChanged();
}

////////////////////////////////////////////////////////
// Routines to help read FITS tables into FTables
////////////////////////////////////////////////////////

TableData* 
FitsTable::loadData(long rowStart, long rowEnd, 
		    const std::vector<std::string>& colMatches) const { 
    // FITS wants column numbers, not names, so try them all one by one
  vector<string> allNames = listFitsColumns(); // index number is column number
  vector<string> colNames;
  vector<int> colNums;
  for (int i=0; i<allNames.size(); i++) {
    set<string> matchout = stringstuff::findMatches(colMatches,
						    vector<string>(1,allNames[i]));
    if (!matchout.empty()) {
      colNames.push_back(allNames[i]);
      colNums.push_back(i);
    }
  }

  // Select row range
  long nrows;
  int status = moveTo();
  fits_get_num_rows(fptr(), &nrows, &status);
  if (rowEnd < 0 ) rowEnd = nrows;
  if (rowEnd < rowStart) rowEnd = rowStart;  // Makes an empty table
  int outRows = rowEnd - rowStart;
  img::TableData* tdata = new TableData(outRows);
  createColumns(tdata, colNames, colNums);
  readFitsData(tdata, colNames, colNums, rowStart, rowEnd);
  return tdata;
}

template <class T>
void addEmptyColumn(img::TableData* tptr, string name, long repeat, long stringLength=-1) {
  if (repeat<0 || repeat>1) tptr->addColumn(vector<vector<T> >(), name, repeat, stringLength);
  else tptr->addColumn(vector<T>(), name, repeat, stringLength);
}

std::vector<std::string>
FitsTable::listFitsColumns() const {
  char colName[FITS::MAX_COLUMN_NAME_LENGTH];

  int status = moveTo();
  int ncols;
  status = fits_get_num_cols(fptr(), &ncols, &status);
  std::vector<std::string> result(ncols);
  do {
    int colNum;
    const char* wildcard="*";
    // Note const_cast to match non-const template for CFITSIO
    fits_get_colname(fptr(), CASEINSEN, 
		     const_cast<char *> (wildcard), colName, &colNum, &status);
    if (status==0 || status==COL_NOT_UNIQUE) {
      Assert(colNum > 0 && colNum <= ncols);
      result[colNum-1] = colName;
    }
  } while (status==COL_NOT_UNIQUE);
  if (status==COL_NOT_FOUND) {
    flushFitsErrors(status);
  }
  checkCFITSIO(status,"findFitsColumns()");
  return result;
}

void
FitsTable::createColumns(img::TableData* tptr, 
			 const vector<string>& colNames,
			 const vector<int>& colNumbers) const {
  int status=0;
  int datatype;
  long repeat, width;
  for (int i=0; i<colNames.size(); i++) {
    status = moveTo();
    fits_get_eqcoltype(fptr(), colNumbers[i]+1, &datatype, &repeat, &width, &status);
    checkCFITSIO(status,"Getting info for column " + colNames[i]);

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
      addEmptyColumn<bool>(tptr, colNames[i], repeat);
      break;
    case FITS::Tlogical:
      addEmptyColumn<bool>(tptr, colNames[i], repeat);
      break;
    case FITS::Tbyte:
      addEmptyColumn<unsigned char>(tptr, colNames[i], repeat);
      break;
    case FITS::Tsbyte:
      addEmptyColumn<signed char>(tptr, colNames[i], repeat);
      break;
    case FITS::Tshort:
      addEmptyColumn<short>(tptr, colNames[i], repeat);
      break;
    case FITS::Tushort:
      addEmptyColumn<unsigned short>(tptr, colNames[i], repeat);
      break;
    case FITS::Tint32bit:
      if (FITS::CIntIsFITSLong)
	// int is 4 bytes, that's what we want
	addEmptyColumn<int>(tptr, colNames[i], repeat);
      else if (!FITS::CLongIsFITSLongLong)
	// If long is not 8 bytes, it's 4 and we'll use it
	addEmptyColumn<long>(tptr, colNames[i], repeat);
      else
	throw FITSError("No intrinsic type holds 32-bit integers!");
      break;
      // Does not appear that unsigned int can be in FITS Table, but put it here in case:
    case FITS::Tuint:
      addEmptyColumn<unsigned int>(tptr, colNames[i], repeat);
      break;
    case FITS::Tulong:
      // ?? ambiguous what this would mean if it arose, go with intrinsic long
      addEmptyColumn<unsigned long>(tptr, colNames[i], repeat);
      break;
    case FITS::Tlonglong:
      addEmptyColumn<LONGLONG>(tptr, colNames[i], repeat);
      break;
    case FITS::Tfloat:
      addEmptyColumn<float>(tptr, colNames[i], repeat);
      break;
    case FITS::Tdouble:
      addEmptyColumn<double>(tptr, colNames[i], repeat);
      break;
    case FITS::Tcomplex:
      addEmptyColumn<complex<float> >(tptr, colNames[i], repeat);
      break;
    case FITS::Tdblcomplex:
      addEmptyColumn<complex<double> >(tptr, colNames[i], repeat);
      break;
    case FITS::Tstring:
      // Note that CFITSIO convention is that an array of fixed-length strings
      // can be specifed, in which case each string is length =width and repeat
      // gives total number of chars in whole array
      // If this was a variable-length char column, then I will call this a 
      // column with a single variable-length string in it:
      if (isVariable) addEmptyColumn<string>(tptr, colNames[i], 1, -1);
      else if (width==repeat) {
	// Each cell is a single fixed-length string
	Assert(repeat>=1);
	addEmptyColumn<string>(tptr, colNames[i], 1, repeat);
      } else {
	// Each cell is a fixed-size array of fixed-length strings
	int nStrings = repeat / width;
	Assert(nStrings>1);
	addEmptyColumn<string>(tptr, colNames[i], nStrings, width);
      }
      break;
    default:
      FormatAndThrow<FITSError>() << "Cannot convert FITS column data type " 
				  << dt << " to intrinsic type";
    } // end switch (dt)
  } // end loop over columns
}

void 
FitsTable::readFitsData(img::TableData* tptr, 
			const vector<string>& names,
			const vector<int>& numbers,
			long rowStart,
			long rowEnd) const {

  // Now read in the data from FITS file into our columns, using the
  // recommended row increments for buffering
  long bufferRows;
  int status = moveTo();
  fits_get_rowsize(fptr(), &bufferRows, &status);
  checkCFITSIO(status, "readFitsData()");
  
  for ( ; rowStart < rowEnd; rowStart+=bufferRows) {
    long rowCount = std::min(rowEnd, rowStart+bufferRows) - rowStart;
    for (int i=0; i<numbers.size(); i++) {
      string name = names[i];
      int num = numbers[i];

      switch ( (*tptr)[name]->elementType() ) {
      case FITS::Tstring:
	getFitsColumnData<string>(tptr, name, num, rowStart, rowCount);
	break;
      case FITS::Tlogical:
	getFitsColumnData<bool>(tptr, name, num, rowStart, rowCount);
	break;
      case FITS::Tbyte:
	getFitsColumnData<unsigned char>(tptr, name, num, rowStart, rowCount);
	break;
      case FITS::Tsbyte:
	getFitsColumnData<signed char>(tptr, name, num, rowStart, rowCount);
	break;
      case FITS::Tshort:
	getFitsColumnData<short>(tptr, name, num, rowStart, rowCount);
	break;
      case FITS::Tushort:
	getFitsColumnData<unsigned short>(tptr, name, num, rowStart, rowCount);
	break;
      case FITS::Tint:
	getFitsColumnData<int>(tptr, name, num, rowStart, rowCount);
	break;
      case FITS::Tuint:
	getFitsColumnData<unsigned int>(tptr, name, num, rowStart, rowCount);
	break;
      case FITS::Tlong:
	getFitsColumnData<long>(tptr, name, num, rowStart, rowCount);
	break;
      case FITS::Tulong:
	getFitsColumnData<unsigned long>(tptr, name, num, rowStart, rowCount);
	break;
      case FITS::Tlonglong:
	getFitsColumnData<LONGLONG>(tptr, name, num, rowStart, rowCount);
	break;
      case FITS::Tfloat:
	getFitsColumnData<float>(tptr, name, num, rowStart, rowCount);
	break;
      case FITS::Tdouble:
	getFitsColumnData<double>(tptr, name, num, rowStart, rowCount);
	break;
      case FITS::Tcomplex:
	getFitsColumnData<complex<float> >(tptr, name, num, rowStart, rowCount);
	break;
      case FITS::Tdblcomplex:
	getFitsColumnData<complex<double> >(tptr, name, num, rowStart, rowCount);
	break;
      default:
	FormatAndThrow<FITSError>() << "Cannot convert data type " 
				    << (*tptr)[name]->elementType() << " to FTable for column "
				    << name;
      }  // end switch on data type

    } // End column loop

  } // End row loop

} // end readFitsData()

// Now a template for routine that gets our column data from file via CFITSIO routines
// Note that this is being called with template being the type of the data in the
// C++ FTable.
template <class T>
void 
FitsTable::getFitsColumnData(img::TableData* tptr, 
			     string colName, int icol, 
			     long rowStart, long nRows) const {
  if (nRows < 1) return;
  FITS::DataType dtype = FITS::FITSTypeOf<T>();
  int status=0;
  T* nullPtr = 0;
  int anyNulls;
  long repeat = (*tptr)[colName]->repeat();
  if (repeat < 0) {
    // Branch for variable-length array: read row by row
    long nElements;
    long offset;
    vector<T> data;
    for (int i=0; i<nRows; i++) {
      status = moveTo();
      // Note adding 1 to row numbers when FITS sees them, since CFITSIO is 1-indexed:
      fits_read_descript(fptr(), icol+1, rowStart+1+i, &nElements, &offset, &status);
      data.resize(nElements);
      fits_read_col(fptr(), dtype, icol+1, (LONGLONG) rowStart+1+i, (LONGLONG) 1, 
		    (LONGLONG) nElements, nullPtr, &data[0], &anyNulls, &status);

      checkCFITSIO(status, "Reading table column " + colName);
      tptr->writeCell(data, colName, rowStart+i);
    } // end row loop
  } else if (repeat>1) {
    // Branch for fixed-length array - can read all rows at one time
    long nElements = repeat * nRows;
    vector<T> data(nElements);
    status = moveTo();
    fits_read_col(fptr(), dtype, icol+1, 
		  (LONGLONG) rowStart+1, (LONGLONG) 1, (LONGLONG) nElements,
		  nullPtr, &data[0], &anyNulls, &status);
    checkCFITSIO(status, "Reading table vector column " + colName);
    typename vector<T>::const_iterator inptr = data.begin();
    for (int i=0; i<nRows; i++) {
      tptr->writeCell(vector<T>(inptr, inptr+repeat), colName, rowStart+i);
      inptr += repeat;
    }
  } else {
    // branch for scalar column
    vector<T> data(nRows);
    status = moveTo();
    fits_read_col(fptr(), dtype, icol+1, 
		  (LONGLONG) rowStart+1, (LONGLONG) 1, (LONGLONG) nRows,
		  nullPtr, &data[0], &anyNulls, &status);
    checkCFITSIO(status, "Reading table column " + colName);
    tptr->writeCells(data, colName, rowStart);
  }
}


// Specialize for logical: CFITSIO is expecting char arrays to hold data, so 
// just do this explicitly instead of assuming that bool = char on this system.
template <>
void 
FitsTable::getFitsColumnData<bool>(img::TableData* tptr, 
				   string colName, int icol, 
				   long rowStart, long nRows) const {
  if (nRows < 1) return;
  FITS::DataType dtype = FITS::Tlogical;
  int status=0;
  char* nullPtr = 0;
  int anyNulls;
  long repeat = (*tptr)[colName]->repeat();
  if (repeat < 0) {
    // Branch for variable-length array: read row by row
    long nElements;
    long offset;
    vector<char> data;
    for (int i=0; i<nRows; i++) {
      status = moveTo();
      // Note adding 1 to row numbers when FITS sees them, since CFITSIO is 1-indexed:
      fits_read_descript(fptr(), icol+1, rowStart+1+i, &nElements, &offset, &status);
      data.resize(nElements);
      fits_read_col(fptr(), dtype, icol+1, (LONGLONG) rowStart+1+i, (LONGLONG) 1, 
		    (LONGLONG) nElements, nullPtr, &data[0], &anyNulls, &status);
      checkCFITSIO(status, "Reading table bool column " + colName);
      tptr->writeCell(vector<bool>(data.begin(), data.end()), colName, rowStart+i);
    } // end row loop
  } else if (repeat>1) {
    // Branch for fixed-length array
    long nElements=repeat * nRows;
    vector<char> data(nElements);
    status = moveTo();
    fits_read_col(fptr(), dtype, icol+1, 
		  (LONGLONG) rowStart+1, (LONGLONG) 1, (LONGLONG) nElements,
		  nullPtr, &data[0], &anyNulls, &status);
    checkCFITSIO(status, "Reading bool table vector column " + colName);
    typename vector<char>::const_iterator inptr = data.begin();
    for (int i=0; i<nRows; i++) {
      tptr->writeCell(vector<bool>(inptr, inptr+repeat), colName, rowStart+i);
      inptr += repeat;
    }
  } else {
    // branch for scalar column
    vector<char> data(nRows);
    status = moveTo();
    fits_read_col(fptr(), dtype, icol+1, 
		  (LONGLONG) rowStart+1, (LONGLONG) 1, (LONGLONG) nRows,
		  nullPtr, &data[0], &anyNulls, &status);
    checkCFITSIO(status, "Reading bool table column " + colName);
    // Convert to a bool vector
    vector<bool> booldata(data.begin(), data.end());
    tptr->writeCells(booldata, colName, rowStart);
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
FitsTable::getFitsColumnData<string>(img::TableData* tptr, string colName, int icol, 
				     long rowStart, long nRows) const {
  if (nRows < 1) return;
  FITS::DataType dtype = Tstring;
  int status=0;
  char* nullPtr = 0;
  int anyNulls;
  long repeat = (*tptr)[colName]->repeat();
  // First consider variable-length strings:
  long length = (*tptr)[colName]->stringLength();
  if (length < 0) {
    if (repeat == 1) {
      // Read a single variable-length string for each row
      long nElements;
      long offset;
      char *data;
      for (int i=0; i<nRows; i++) {
	status = moveTo();
	// Note adding 1 to row numbers when FITS sees them, since CFITSIO is 1-indexed:
	fits_read_descript(fptr(), icol+1, rowStart+1+i, &nElements, &offset, &status);
	data = new char[nElements+1];
	// Read as a c-string; wants char** for destination
	fits_read_col(fptr(), dtype, icol+1, 
		      (LONGLONG) rowStart+1+i, (LONGLONG) 1, (LONGLONG) nElements, 
		      nullPtr, &data, &anyNulls, &status);
	checkCFITSIO(status, "Reading string table column " + colName);
	// Rewrite as string: make sure it's null-terminated
	data[nElements] = 0;
	string s(data);
	tptr->writeCell(s, colName, rowStart+i);
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
    status = moveTo();
    // ??? nStrings or # characters here ???
    fits_read_col(fptr(), dtype, icol+1, 
		  (LONGLONG) rowStart+1, (LONGLONG) 1, (LONGLONG) nStrings,
		  nullPtr, &data[0], &anyNulls, &status);
    checkCFITSIO(status, "Reading table column " + colName);

    if (repeat==1) {
      // One string per cell; make an array of them all
      vector<string> vs(nRows);
      for (int i=0; i<nRows; i++) {
	// Make sure it's null terminated:
	data[i][length]='\0';
	vs[i] = data[i];
	delete[] data[i];
      }
      tptr->writeCells(vs, colName, rowStart);
    } else {
      // Multiple strings per cell: write cells one at a time
      vector<string> vs(repeat);
      int iString = 0;
      for (int j=0; j<nRows; j++) {
	for (int i=0; i<repeat; i++, iString++) {
	  data[iString][length] = '\0'; // insure null termination
	  vs[i] = data[iString];
	  delete[] data[iString];
	}
	tptr->writeCell(vs, colName, rowStart+j);
      }
    }
  } // end if/then/else for array configuration
}


////////////////////////////////////////////////////////////////////////////////
////////// Writing routines
////////////////////////////////////////////////////////////////////////////////

// Write contents of FTable into the FitsTable (not header)
void
FitsTable::writeFitsTableData(const img::TableData *tptr) {
  clearFitsTableData();
  long rowStart = 0;
  long rowEnd = tptr->nrows();
  vector<string> colNames;
  vector<int> colNums;
  addFitsColumns(tptr, colNames, colNums);
  // Now going to write out the contents to FITS, using recommended buffering row intervals.
  long bufferRows;
  int status = moveTo();
  fits_get_rowsize(fptr(), &bufferRows, &status);
  checkCFITSIO(status, "writeFitsTableData()"); 
  
  for ( ; rowStart < rowEnd; rowStart+=bufferRows) {
    long rowCount = std::min(rowEnd, rowStart+bufferRows) - rowStart;
    // Loop over columns, writing interval of each to output
    for (int iCol = 0; iCol < colNames.size(); iCol++) {
      FITS::DataType dt = (*tptr)[colNames[iCol]]->elementType();
      // Dispatch the data according to the datatype:
      switch (dt) {
      case Tlogical:
	writeFitsColumn<bool>(tptr, colNames[iCol], colNums[iCol], 
			      rowStart, rowStart+rowCount);
	break;
      case Tstring:
	writeFitsColumn<string>(tptr, colNames[iCol], colNums[iCol], 
				rowStart, rowStart+rowCount);
	break;
      case Tbyte:
	writeFitsColumn<unsigned char>(tptr, colNames[iCol], colNums[iCol], 
				       rowStart, rowStart+rowCount);
	break;
      case Tsbyte:
	writeFitsColumn<signed char>(tptr, colNames[iCol], colNums[iCol], 
				     rowStart, rowStart+rowCount);
	break;
      case Tshort:
	writeFitsColumn<short>(tptr, colNames[iCol], colNums[iCol], 
			       rowStart, rowStart+rowCount);
	break;
      case Tushort:
	writeFitsColumn<unsigned short>(tptr, colNames[iCol], colNums[iCol], 
					rowStart, rowStart+rowCount);
	break;
      case Tint:
	writeFitsColumn<int>(tptr, colNames[iCol], colNums[iCol], 
			     rowStart, rowStart+rowCount);
	break;
      case Tuint:
	writeFitsColumn<unsigned int>(tptr, colNames[iCol], colNums[iCol], 
				      rowStart, rowStart+rowCount);
	break;
      case Tlong:
	writeFitsColumn<long>(tptr, colNames[iCol], colNums[iCol], 
			     rowStart, rowStart+rowCount);
	break;
      case Tulong:
	writeFitsColumn<unsigned long>(tptr, colNames[iCol], colNums[iCol], 
				      rowStart, rowStart+rowCount);
	break;
      case Tlonglong:
	writeFitsColumn<LONGLONG>(tptr, colNames[iCol], colNums[iCol], 
				  rowStart, rowStart+rowCount);
	break;
      case Tfloat:
	writeFitsColumn<float>(tptr, colNames[iCol], colNums[iCol], 
			       rowStart, rowStart+rowCount);
	break;
      case Tdouble:
	writeFitsColumn<double>(tptr, colNames[iCol], colNums[iCol], 
				rowStart, rowStart+rowCount);
	break;
      case Tcomplex:
	writeFitsColumn<complex<float> >(tptr, colNames[iCol], colNums[iCol], 
					 rowStart, rowStart+rowCount);
	break;
      case Tdblcomplex:
	writeFitsColumn<complex<double> >(tptr, colNames[iCol], colNums[iCol], 
					  rowStart, rowStart+rowCount);
	break;
      default:
	ostringstream oss;
	oss << "Unknown DataType " << dt
	    << " for writing data to FITS column";
	throwFitsOrDump(oss.str());
      } // end DataType switch
    } // end column loop
  } // end row chunk loop
}

void
FitsTable::addFitsColumns(const img::TableData* tptr, 
			  vector<string>& colNames, vector<int>& colNums) {
  checkWriteable("addFitsColumns()");
  int colNum;
  // Add these columns to end of any there:
  int status = moveTo();
  fits_get_num_cols(fptr(), &colNum, &status);
  checkCFITSIO(status, "addFitsColumns()");

  colNums.clear();

  // Get all column names
  colNames.clear();
  for (img::TableData::const_iterator i = tptr->begin();
       i != tptr->end();
       ++i) 
    colNames.push_back((*i)->name());

  for (int i=0; i<colNames.size(); i++, colNum++) {
    // loop over ft's columns:
    // First make the TFORM string
    long repeat = (*tptr)[colNames[i]]->repeat();
    ostringstream oss;
    if ( (*tptr)[colNames[i]]->elementType() == FITS::Tstring) {
      long length = (*tptr)[colNames[i]]->stringLength();
      if (length < 0) {
	// single variable-length string per cell:
	if (repeat==1) oss << "1PA";
	// otherwise: we will punt right now
	else
	  throwFitsOrDump("Cannot format array of variable-length strings for bintable "
			  "at column " + colNames[i]);
      } else {
	// Fixed-length strings
	if (repeat < 0) {
	  throwFitsOrDump("Cannot format variable-length arrays of strings "
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
      oss << (*tptr)[colNames[i]]->columnCode();
    }

    status = moveTo();
    fits_insert_col(fptr(), colNum+1, const_cast<char*> (colNames[i].c_str()), 
		    const_cast<char*> (oss.str().c_str()),
		    &status);
    colNums.push_back(colNum);
    checkCFITSIO(status, "creating column " + colNames[i]); 
  }
}

void
FitsTable::clearFitsTableData() {
  // Remove all data from Fits Table extension
  checkWriteable("in clearFitsTableData()");

  int status = moveTo();
  int ncols;
  // Delete all columns, from back to front:
  fits_get_num_cols(fptr(), &ncols, &status);
  for (int i=ncols; i>0; i--)
    fits_delete_col(fptr(), i, &status);
  // compress (=empty) the heap:
  fits_compress_heap(fptr(), &status);
  checkCFITSIO(status,"In FitsTable::clear(): ");
}

template <class DT> 
void
FitsTable::writeFitsColumn(const img::TableData *tptr, 
			   string colName, int colNum,
			   long rowStart, long rowEnd) {
  checkWriteable("writeFitsColumn()");
  int status=0;
  long repeat = (*tptr)[colName]->repeat();
  vector<DT> data;
  if (repeat==1) {
    // Scalar column, can write whole group at once
    tptr->readCells(data, colName, rowStart, rowEnd);
    status = moveTo();
    fits_write_col(fptr(), FITS::FITSTypeOf<DT>(), colNum+1, (LONGLONG) rowStart+1,
		   (LONGLONG) 1, (LONGLONG) (rowEnd-rowStart), &data[0], &status);
    checkCFITSIO(status, "Writing to FITS column " + colName);
  } else {
    // Array cells, write 1 row at a time
    for (int iRow = rowStart; iRow < rowEnd; iRow++) {
      tptr->readCell(data, colName, iRow);
      if (repeat >=0) //**Assert(data.size()==repeat);
	/**/if (data.size() != repeat) 
	  /**/	  cerr << "Writing column " << colName << " to extension " << getName()
	       << " of " << getFilename() 
	       << " row " << iRow 
	       << " data.size() " << data.size()
	       << " repeat " << repeat
	       << endl;
      status = moveTo();
      fits_write_col(fptr(), FITS::FITSTypeOf<DT>(), colNum+1, (LONGLONG) iRow+1,
		     (LONGLONG) 1, (LONGLONG) data.size(), &data[0], &status);
      checkCFITSIO(status, "Writing to FITS column " + colName);
    } // end row loop
  } // End if/else for array/scalar cells
}

// Specialize for logical, need to pass data in/out as char
template <> 
void
FitsTable::writeFitsColumn<bool>(const img::TableData* tptr, 
				 string colName, int colNum,
				 long rowStart, long rowEnd) {
  checkWriteable("Attempt to write to read-only HDU in FITS file " + getFilename());
  int status=0;
  long repeat = (*tptr)[colName]->repeat();
  vector<bool> data;
  if (repeat==1) {
    // Scalar column, can write whole group at once
    tptr->readCells(data, colName, rowStart, rowEnd);
    char cdata[data.size()];
    for (int i=0; i<data.size(); i++)
      cdata[i] = data[i] ? 1 : 0;
    status = moveTo();
    fits_write_col(fptr(), FITS::Tlogical, colNum+1, (LONGLONG) rowStart+1,
		   (LONGLONG) 1, (LONGLONG) (rowEnd-rowStart), cdata, &status);
    checkCFITSIO(status, "Writing to FITS column " + colName);
  } else {
    // Array cells, write 1 row at a time
    for (int iRow = rowStart; iRow < rowEnd; iRow++) {
      tptr->readCell(data, colName, iRow);
      char cdata[data.size()];
      for (int i=0; i<data.size(); i++)
	cdata[i] = data[i] ? 1 : 0;
      if (repeat >=0) Assert(data.size()==repeat);
      status = moveTo();
      fits_write_col(fptr(), FITS::Tlogical, colNum+1, (LONGLONG) iRow+1,
		     (LONGLONG) 1, (LONGLONG) data.size(), cdata, &status);
      checkCFITSIO(status, "Writing to FITS column " + colName);
    } // end row loop
  } // End if/else for array/scalar cells
}

// Specialize for strings
template<>
void
FitsTable::writeFitsColumn<string>(const img::TableData* tptr, 
				   string colName, int colNum,
				   long rowStart, long rowEnd) {
  checkWriteable("writeFitsColumn<string>()");
  int status=0;
  long repeat = (*tptr)[colName]->repeat();
  long length = (*tptr)[colName]->stringLength();
  if (repeat==1) {
    // One string per cell.  Read individually, whether or not fixed-length: FITSIO will deal
    string data;
    for (int iRow = rowStart; iRow < rowEnd; iRow++) {
      tptr->readCell(data, colName, iRow);
      char* cdata[1];
      cdata[0] = const_cast<char*> (data.c_str());
      status = moveTo();
      // fits_write_col wants char** input for source of data
      fits_write_col(fptr(), Tstring, colNum+1, (LONGLONG) iRow+1,
		     (LONGLONG) 1, (LONGLONG) 1,
		     cdata, &status);
      checkCFITSIO(status, "Writing to FITS column " + colName);
    }
  } else {
      if (length < 0 || repeat < 0)
	throwFitsOrDump("Cannot write column <" + colName + "> to FITS because"
			" it has variable-length string arrays.");
      // Array cells, write 1 row at a time
      vector<string> data;
      char* cdata[repeat];
      for (int iRow = rowStart; iRow < rowEnd; iRow++) {
	tptr->readCell(data, colName, iRow);
	data.resize(repeat);
	for (int j=0; j<repeat; j++)
	  cdata[j] = const_cast<char*> (data[j].c_str());
	status = moveTo();
	fits_write_col(fptr(), Tstring, colNum+1, (LONGLONG) iRow+1,
		       (LONGLONG) 1, (LONGLONG) repeat, cdata, &status);
	checkCFITSIO(status, "Writing to FITS column " + colName);
    } // end row loop
  } // End if/else for array/scalar cells
}



