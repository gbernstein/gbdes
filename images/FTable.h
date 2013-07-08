// Analogize the Image class to a 2d table that is an in-memory representation of a FITS
// ASCII or binary table.  Will use the same header structure as the Image class.

// The class is called FTable to distinguish from the Table class used for interpolation.
// The corresponding FITS file HDU's will be called FITSTable class.

// Note that convention here is that rows are 0-indexed, like C++ containers,
// but FITS/CFITSIO has convention of 1-indexed rows.
// Also that insertion / erase indexing follows the C++ container conventions.

#ifndef FTABLE_H
#define FTABLE_H

#include <vector>
#include <set>

#include "Header.h"

namespace img {

  // Exception classes:
  class FTableError: public std::runtime_error {
  public: 
    FTableError(const string &m=""): std::runtime_error("FTable Error: " + m) {}

  };
  class FTableNonExistentColumn: public FTableError {
  public: 
    FTableNonExistentColumn(const string colname=""): 
      FTableError("Non-existent column name <" + colname + ">") {}
  };
  class FTableRangeError: public FTableError {
  public: 
    FTableRangeError(const string& m=""): 
      FTableError("Row out of range " + m) {}
  };
  class FTableLockError: public FTableError {
  public: 
    FTableLockError(const string& m=""): 
      FTableError("Write access to locked data " + m) {}
  };

}

// Include classes that the user will not need to use directly:
#include "FTable2.h"

namespace img {
  // This is a handle class that contains a (reference-counted) header and table
  // Locking occurs at the header and data classes and will be patrolled there.
  class FTable {
  private:
    TableData* D;	//pixel data
    mutable int* dcount;  // link count for the data structure
    Header* H;	//the "header"
    mutable int* hcount;

  public:
    // Constructor for new table that has no data, reserve nrow spaces in all new columns
    FTable(long rowsReserved=0):
      D(new TableData(rowsReserved)),
      H(new Header()),
      dcount(new int(1)),
      hcount(new int(1)) {}
    // Copy constructor is really a data share; lock status is preserved
    // by sharing so this can become unwriteable if we op= it to a locked rhs
    FTable(const FTable &rhs): 
      D(rhs.D),
      H(rhs.H), 
      dcount(rhs.dcount),
      hcount(rhs.hcount) {(*dcount)++; (*hcount)++;}
    // So is assignment:
    FTable& operator=(const FTable& rhs) {
      if (&rhs == this) return *this;
      if (D!=rhs.D) {
	if (--(*dcount)==0) {delete D; delete dcount;}
	D = rhs.D; dcount=rhs.dcount; (*dcount)++;
	D->touch(); // ??? can omit this ???
      }
      if (H!=rhs.H) {
	if (--(*hcount)==0) {delete H; delete hcount;}
	H = rhs.H; hcount=rhs.hcount; (*hcount)++;
	H->touch();  // ??? can omit this ???
      }
      return *this;
    }

    // Constructor that only FitsTable should use, build straight from elements:
    FTable(Header* hh, int* hc, TableData* dd, int *dc):
      D(dd), dcount(dc), H(hh), hcount(hc) {
      (*dcount)++;
      (*hcount)++;
    }
    // Also expose guts: should only use this in FitsTable
    void showGuts(Header* &hh, int* &hc, TableData* &dd, int* &dc) const {
      hh = H; hc = hcount; dd = D; dc = dcount;
    }

    // Create a fresh deep copy
    FTable duplicate() const {
      FTable dup;  // Make empty header & table data
      *(dup.H) = *H;
      delete dup.D;  // clean out its empty TableData
      dup.D = D->duplicate(); // and replace with duplicate of this->D.
      return dup;
    }

    ~FTable() {
      if (--(*dcount)==0) {delete D; delete dcount;}
      if (--(*hcount)==0) {delete H; delete hcount;}
    }

    //////////////////////////////////////
    // Table / column information
    //////////////////////////////////////

    bool isChanged() const {return H->isChanged() || D->isChanged();}
    void clearChanged() {H->clearChanged(); D->clearChanged();}
    bool isLocked() const {return H->isLocked() || D->isLocked();}
    // Note ***there is no unlocking of data.***  This avoids becoming unconst.
    void setLock() {H->setLock(); D->setLock();}


    // Access header:
    Header* header() {return H;}
    const Header* header() const {return H;}

    // Report current table rows
    long nrows() const {return D->nrows();}

    // Report number of columns
    int ncols() const {return D->ncols();}

    // Report names of columns
    vector<string> listColumns() const {return D->listColumns();}

    // Get repeat count of a column's cell: 
    // 1 for scalar, -1 for variable-length, >=0 for fixed.
    long repeat(string columnName) const {return (*D)[columnName]->repeat();}
   
    // If the column holds strings, return the defined max length of strings, 
    // or -1 for variable-length strings.
    long stringLength(string columnName) const {
      return (*D)[columnName]->stringLength();
    }

    // Get the type of data stored in the C arrays (using the FitsTypes.h class)
    FITS::DataType elementType(string columnName) const {
      return (*D)[columnName]->elementType();
    }

    //////////////////////////////////////
    // Reads and writes
    //////////////////////////////////////

    void clear() {H->clear(); D->clear();}

    // Access single element (copy created)
    template <class T>
    void readCell(T& value, string columnName, long row) const {
      D->readCell(value, columnName, row);
    }
    // Or a range: rowEnd is one-past-last row; -1 means go to end.
    template <class T>
    void readCells(vector<T>& values, string columnName, 
		   long rowStart=0, long rowEnd=-1) const {
      D->readCells(values, columnName, rowStart, rowEnd);
    }
    
    // Set elements.  Table grows (with default-initialized data) if row is past current end.
    // When setting fixed-length array cells, too-short array is padded with default value.
    // Too-long input array throws exception.  Same thing for strings, which are treated 
    // as scalar cells in this code but are arrays in FITS BinTable.
    template <class T>
    void writeCell(const T& value, string columnName, long row) {
      D->writeCell(value, columnName, row);
    }
    // Or a range: length of input vector determines rows altered; grow table if needed.
    template <class T>
    void writeCells(const vector<T>& values, string columnName, long rowStart=0) {
      D->writeCells(values, columnName, rowStart);
    }


    //////////////////////////////////////
    // Row and column manipulations
    //////////////////////////////////////

    // Add column from vector of values.  Grow table if needed, pad input array with default
    // value if it is too short.
    // The repeat argument has meaning for string or array columns: <0 for variable length,
    // >=0 gives fixed length.
    // Given stringLength>=0 means all strings must be at this length or shorter (and
    // they will be stored as fixed-length in FITS files).
    // Exception thrown input arrays/strings are too long for fixed-length fields
    template<class T>
      void addColumn(const vector<T>& values, string columnName, long repeat=-1,
		     long stringLength=-1 ) {
      D->addColumn(values, columnName, repeat, stringLength);
    }

    // ??? get/set the column attributes: units, TDIM, format, null values?

    // Delete column: exception to delete non-existent column
    void eraseColumn(string columnName) {D->erase(columnName);}

    // Delete row(s) - not particularly efficient since using vector storage
    // Asking to erase past end is ignored; erase before beginning throws exception.
    void eraseRow(long row) {D->eraseRows(row, row+1);}
    //  rowEnd=one-past-end; -1= go to end.  
    void eraseRows(long rowStart, long rowEnd=-1) {D->eraseRows(rowStart, rowEnd);}

    // Insert rows filled with default values in every column.  Not very efficient.
    // insertBeforeRow = -1 will add row(s) to the end.  Exception to insert beyond end.
    void insertRows(long insertBeforeRow, long insertNumber) {
      D->insertRows(insertBeforeRow, insertNumber);
    }

    // Get new table with rows chosen various ways:
    // range of row numbers; also keep only columns that match one of the regexp's:
    FTable extract(long rowStart=0, long rowEnd=-1,
		   const vector<string>& regexps
		   = vector<string>(1,".*")) const;

    // bool that signals inclusion of a row
    FTable extractRows(vector<bool>& vb) const;

    // expression to evaluate for a row
    FTable extractRows(const string& expression) const;

    // rows with column's value in container
    // column data should be convertible to Container::key_type 
    template <class Container>
    FTable extractRows(const string& columnName, const Container& keepers) const {
      FTable result;
      D->filterRows(result.D, columnName, keepers);
      result.D->clearChanged();
      result.H->copyFrom(*H);
      return result;
    } 

    // Same as above but change internally rather than issuing new FTable
    void filter(long rowStart=0, long rowEnd=-1,
		const vector<string>& regexps
		= vector<string>(1,".*"));
    // bool that signals inclusion of a row
    void filterRows(vector<bool>& vb);
    // expression to evaluate for a row
    void filterRows(const string& expression);
    // rows with column's value in container
    template <class Container>
    void filterRows(const string& columnName, const Container& keepers) {
      D->filterRows(D, columnName, keepers);
    }

    //////////////////////////////////////
    // Evaluate arithmetic / logical expression on columns, 
    // convert to type T and place in vector
    template <class T>
    void evaluate(vector<T>& result, string expression) {
      D->evaluate(result, expression, H);
    }

    //////////////////////////////////////
    // Sort in increasing order with native comparison for specified column:
    void sort(string columnName);
    // Or according to specified ordering function operating on that column:
    template <class LessThan>
    void sort(string columnName, LessThan& ordering);

    ///////////////////////
    // Direct header access, for convenience:
    // Get/set the value of header records.  Bool returns false if
    // keyword doesn't exist or does not match type of argument.
    template <class U> 
    bool getHdrValue(const string keyword, U& outVal) const {
      return header()->getValue(keyword, outVal);
    }
    template <class U> 
    bool setHdrValue(const string keyword, const U& inVal) {
      return header()->setValue(keyword,inVal);
    }

  };

} // namespace img

#endif // FTABLE_H
