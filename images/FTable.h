// $Id: FTable.h,v 1.2 2012/09/08 15:48:19 garyb Exp $ 
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

#include "Image.h"


namespace img {


  // Exception classes:
  class FTableError: public MyException {
  public: 
    FTableError(const string &m=""): MyException("FTable Error: " + m) {}

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

}

// Include classes that the user will not need to use directly:
#include "FTable2.h"

namespace img {
  // This is a handle class that contains a (reference-counted) header and table
  class FTable {
  private:
    TableData* D;	//pixel data
    mutable int* dcount;  // link count for the data structure
    ImageHeader* H;	//the "header"
    mutable int* hcount;

  public:
    // Constructor for new table that has no data, reserve nrow spaces in all new columns
    FTable(long rowsReserved=0):
      D(new TableData(rowsReserved)),
      H(new ImageHeader()),
      dcount(new int(1)),
      hcount(new int(1)) {}
    // Copy constructor is really a data share:
    FTable(const FTable &rhs): // ??? how to keep from setting non-const to const
      D(rhs.D),
      H(rhs.H), 
      dcount(rhs.dcount),
      hcount(rhs.hcount) {(*dcount)++; (*hcount)++;}
    // So is assignment:
    FTable& operator=(const FTable& rhs) {
      // Note no assignment of const image to non-const image. ???
      if (&rhs == this) return *this;
      if (D!=rhs.D) {
	if (--(*dcount)==0) {delete D; delete dcount;}
	D = rhs.D; dcount=rhs.dcount; (*dcount)++;
      }
      if (H!=rhs.H) {
	if (--(*hcount)==0) {delete H; delete hcount;}
	H = rhs.H; hcount=rhs.hcount; (*hcount)++;
      }
      return *this;
    }
    // Create a fresh deep copy
    FTable duplicate() const {
      FTable dup;  // Make empty header & table data
      *(dup.H) = *H;
      delete dup.D;  // clean out its empty TableData
      dup.D = D->duplicate(); // and replace with duplicate of this->D.
      return dup;
    }
    // Or adopt another table's header only:
    void copyHeader(const ImageHeader& rhs) {*H=rhs;}

    ~FTable() {
      if (--(*dcount)==0) {delete D; delete dcount;}
      if (--(*hcount)==0) {delete H; delete hcount;}
    }

    // Access elements directly:
    ImageHeader* header() {return H;}
    const ImageHeader* header() const {return H;}
    TableData* data() {return D;}
    const TableData* data() const {return D;}

    // Access single element (copy created)
    template <class T>
    void readCell(T& value, string columnName, long row) const {
      D->readCell(value, columnName, row);
    }
    // Or a range: rowEnd is one-past-last row; -1 means go to end.
    template <class T>
    void readCells(vector<T>& values, string columnName, long rowStart=0, long rowEnd=-1) const {
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

    // Report current table rows
    long nrows() const {return D->nrows();}

    // Report number of columns
    int ncols() const {return D->ncols();}

    // Report names of columns
    vector<string> columnNames() const {
      vector<string> out;
      for (TableData::const_iterator i = D->begin();
	   i != D->end();
	   ++i) 
	out.push_back((*i)->name());
      return out;
    }

    // Get repeat count of a column: 1 (or 0) for scalar, -1 for variable-length, >=0 for fixed.
    long repeat(string columnName) const {return (*D)[columnName]->repeat();}
    // If the column holds strings, return the defined max length of strings, or -1 for variable.
    long stringLength(string columnName) const {return (*D)[columnName]->stringLength();}
    // Get the type of data stored in the C arrays (using the FITStypes.h class)
    FITS::DataType elementType(string columnName) const {return (*D)[columnName]->elementType();}

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
    void deleteColumn(string columnName) {D->erase(columnName);}

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

    // Evaluate arithmetic / logical expression on columns, convert to type T and place in vector
    template <class T>
    void evaluate(vector<T>& result, string expression);
    // Delete all rows that have zero in the input vector
    template <class T>
    void filter(const vector<T>& retain);
    // Or filter those that give true (non-zero) result in expression
    void filter(string expression);

    // Sort in increasing order with native comparison for specified column:
    void sort(string columnName);
    // Or according to specified ordering function operating on that column:
    template <class LessThan>
    void sort(string columnName, LessThan& ordering);

    // Retain only elements whose keyword matches an element of the set:
    template <class T>
    void select(string columnName, const set<T>& retain);

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
