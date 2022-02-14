// In-memory version of the FITS table structure.
// This file contains classes / declarations that should not be needed
// by the user.

// Define RANGE_CHECK to have row arguments always checked:
#define RANGE_CHECK //???

#ifndef FTABLE2_H
#define FTABLE2_H

#include <vector>
#include <map>
#include <typeinfo>

#include "StringStuff.h"

// sorts & filters

namespace img {

  class ColumnBase {
  public:
    ColumnBase(string name): colName(name) {}
    virtual ~ColumnBase() {}
    virtual ColumnBase* duplicate() const=0;
    const string& name() const {return colName;}
    // Number of rows of data
    virtual long size() const =0;
    virtual void reserve(long rowsReserved) =0;
    virtual void resize(long newSize) =0;

    // Return 1 (or 0) for scalar cells, >1 for fixed-length array cell, <0 for var-length
    virtual long repeat() const=0;
    // Return max allowed length of string-valued elements of column, or -1 to be free
    virtual long stringLength() const=0;
    // FITS DataType for scalar or array element:
    virtual FITS::DataType elementType() const =0;  
    // FITSIO code for appropriate column data storage:
    virtual char columnCode() const =0;  

    virtual void eraseRows(long rowStart, long rowEnd) =0;
    // insertBeforeRow = -1 will add row(s) to the end.
    // Return value is number of rows after insertion.
    virtual long insertRows(long insertBeforeRow, long insertNumber) =0;

    // Create new column that is a subset of this one:
    virtual ColumnBase* copyRows(long rowStart, long rowEnd) const =0;
    // Choose the rows with true in their column.  
    // nReserve is size of vector needed, -1 means don't know
    virtual ColumnBase* copyRows(std::vector<bool>& vb,
				 int nReserve=-1) const =0;
    // Convert column data to requested vector, if no loss of information
    // Can't template this because it's a virtual function!
    virtual bool makeVector(std::vector<bool>& vout) const {return false;}
    virtual bool makeVector(std::vector<string>& vout) const {return false;}
    virtual bool makeVector(std::vector<long>& vout) const {return false;}
    virtual bool makeVector(std::vector<double>& vout) const {return false;}
  private:
    const string colName;
    // Ancillary information ???
  };

  template <class DT>
  class ScalarColumn: public ColumnBase {
  protected:
    vector<DT> v;  // Holds the column's data
    // Protected function to check length of input strings
    void checkLength(const string& input) {
      if (input.size() > length) 
	throw FTableError("Input string <" + input + "> too long at column " + name());
    }
  private:
    const FITS::DataType dType;
    const char colCode;
    long length;   // max length of string, if this is holding strings
  public:
    ScalarColumn(string name, long length_=-1): 
      ColumnBase(name), dType(FITS::FITSTypeOf<DT>()), 
      colCode(FITS::ColumnCode<DT>()), length(length_) {}
    ScalarColumn(const vector<DT>& in, string name, long length_=-1):
      ColumnBase(name), dType(FITS::FITSTypeOf<DT>()), 
      colCode(FITS::ColumnCode<DT>()), v(in), length(length_) {}
    virtual ~ScalarColumn() {}
    virtual ColumnBase* duplicate() const {return new ScalarColumn(*this);}
    virtual long size() const {return v.size();}
    virtual void reserve(long rowsReserved) {v.reserve(rowsReserved);}
    virtual void resize(long newSize) {v.resize(newSize);}
    virtual long repeat() const {return 1;}  // scalar default
    virtual long stringLength() const {return length;} 
    virtual FITS::DataType elementType() const {return dType;}
    virtual char columnCode() const {return colCode;}

    // Convert column data to requested vector, if no loss of information
    // Can't template this because it's a virtual function!
    // Will specialize some bool and string
    virtual bool makeVector(std::vector<bool>& vout) const {return false;}
    virtual bool makeVector(std::vector<string>& vout) const {return false;}
    virtual bool makeVector(std::vector<long>& vout) const {return false;}
    virtual bool makeVector(std::vector<double>& vout) const {return false;}

    // Create new column that is a subset of this one:
    virtual ColumnBase* copyRows(long rowStart, long rowEnd) const {
      long nRows = rowEnd - rowStart;
      ScalarColumn<DT>* target = new ScalarColumn<DT>(name(), stringLength());
      if (nRows>0) target->v.resize(nRows);
      for (int i=0; i<nRows; i++)
	target->v[i] = v[i+rowStart];
      return target;
    }
    virtual ColumnBase* copyRows(std::vector<bool>& vb, int nReserve) const {
      if (vb.size() > v.size())
	throw FTableError("too many input bools in copyRows() on column " + name());
      ScalarColumn<DT>* target = new ScalarColumn<DT>(name(), stringLength());
      if (nReserve<0) {
	nReserve = 0;
	for (std::vector<bool>::const_iterator i=vb.begin();
	     i != vb.end();
	     ++i)
	  if (*i) nReserve++;
      }
      target->v.resize(nReserve);
      typename vector<DT>::iterator j = target->v.begin();
      for (int i=0; i<vb.size(); i++)
	if (vb[i]) *(j++) = v[i];
      return target;
    }

    // Range checking will be done at TableData level, so none here
    virtual void readCell(DT& value, long row) const {value = v[row];}
    virtual void readCells(vector<DT>& values, long rowStart=0, long rowEnd=-1) const {
      int nCopy = (rowEnd >=0 ? rowEnd : v.size()) - rowStart;
      if (nCopy < 0) nCopy = 0;
      values.resize(nCopy);
      for (int i=0; i<nCopy; i++)
	values[i] = v[i+rowStart];
    }

    // Writes return the length of column after writing.
    virtual long writeCell(const DT& value, long row) {
      // extend column to hold data
      if (row >= v.size()) v.resize(row+1);
      v[row] = value;
      return v.size();
    }
    // Or a range: length of input vector determines rows altered; grow table if needed.
    virtual long writeCells(const vector<DT>& values, long rowStart=0) {
      long requiredSize = rowStart + values.size();
      if (requiredSize > v.size()) v.resize(requiredSize);
      for (int i=0; i<values.size(); i++)
	v[i+rowStart] = values[i];
      return v.size();
    }

    virtual void eraseRows(long rowStart, long rowEnd) {
      Assert(rowStart > 0);
      Assert(rowEnd > rowStart);
      Assert(rowEnd <= v.size());
      v.erase(v.begin()+rowStart, v.begin()+rowEnd);
    }

    // insertBeforeRow = -1 will add row(s) to the end.
    virtual long insertRows(long insertBeforeRow, long insertNumber) {
      if (insertBeforeRow < 0 || insertBeforeRow >= v.size())
	v.insert(v.end(), insertNumber, DT());
      else
	v.insert(v.begin()+insertBeforeRow, insertNumber, DT());
      return v.size();
    }
    
  };

  // Specialization for string-valued scalar arrays, check length as needed
  template<>
  inline ScalarColumn<string>::ScalarColumn(const vector<string>& in, string name, 
				     long length_): 
    ColumnBase(name), dType(FITS::FITSTypeOf<string>()), 
    colCode(FITS::ColumnCode<string>()), length(length_) {
    if (stringLength()>=0) 
      for (int i=0; i<in.size(); i++)
	checkLength(in[i]);
    v = in;
  }

  template<>
  long
  inline ScalarColumn<string>::writeCell(const string& value, long row) {
    if (stringLength()>=0) checkLength(value);
    // extend column to hold data
    if (row >= v.size()) v.resize(row+1);
    v[row] = value;
    return v.size();
  }

  template<>
  long
  inline ScalarColumn<string>::writeCells(const vector<string>& values, long rowStart) {
    if (stringLength()>=0) 
      for (int i=0; i<values.size(); i++)
	checkLength(values[i]);
    long requiredSize = rowStart + values.size();
    if (requiredSize > v.size()) v.resize(requiredSize);
    for (int i=0; i<values.size(); i++)
      v[i+rowStart] = values[i];
    return v.size();
  }
  
  ////////////////////////
  // Specializations of the vector-making routine
#define MVECTOR(T1,T2)                                         \
   template<>						       \
   inline bool	  				               \
   ScalarColumn<T2>::makeVector(std::vector<T1>& vout) const { \
     vout.resize(v.size());				       \
     for (int i=0; i<v.size(); i++) vout[i]=v[i];              \
     return true;					       \
   }							       

  MVECTOR(bool, bool)
  MVECTOR(string, string)
  MVECTOR(long, char)
  MVECTOR(long, signed char)
  MVECTOR(long, short)
  MVECTOR(long, unsigned short)
  MVECTOR(long, int)
  MVECTOR(long, unsigned int)
  MVECTOR(long, long)
  MVECTOR(long, unsigned long)
  MVECTOR(double, char)
  MVECTOR(double, signed char)
  MVECTOR(double, short)
  MVECTOR(double, unsigned short)
  MVECTOR(double, int)
  MVECTOR(double, unsigned int)
  MVECTOR(double, long)
  MVECTOR(double, unsigned long)
  MVECTOR(double, float)
  MVECTOR(double, double)
#undef MVECTOR

  template <class DT>
  class ArrayColumn: public ScalarColumn<vector<DT> > {
  private:
    typedef ScalarColumn<vector<DT> > Base;
    const FITS::DataType dType2;
    const char colCode2;
  protected:
    using ScalarColumn<vector<DT> >::v;
  public:
    ArrayColumn(string name, long length_=-1): 
      Base(name, length_), dType2(FITS::FITSTypeOf<DT>()),
      colCode2(FITS::ColumnCode<DT>()) {}
    ArrayColumn(const vector<vector<DT> >& in, string name, long length_=-1): 
      Base(in, name, length_), dType2(FITS::FITSTypeOf<DT>()),
      colCode2(FITS::ColumnCode<DT>()) {}
    virtual ColumnBase* duplicate() const {return new ArrayColumn(*this);}
    virtual ~ArrayColumn() {}

    // For variable-length array, the only thing we need to override is 
    // the repeat() count and the elementType should be DT, not the vector<DT>:
    virtual long repeat() const {return -1;}
    virtual FITS::DataType elementType() const {return dType2;}
    virtual char columnCode() const {return colCode2;}

    // Need to declare these so they can be specialized for string:
    virtual long writeCell(const vector<DT>& value, long row) {
      return Base::writeCell(value,row);
    }
    virtual long writeCells(const vector<vector<DT> >& values, long rowStart=0) {
      return Base::writeCells(values, rowStart);
    }

    // Create new column that is a subset of this one:
    virtual ColumnBase* copyRows(long rowStart, long rowEnd) const {
      long nRows = rowEnd - rowStart;
      ArrayColumn<DT>* target = 
	new ArrayColumn<DT>(ColumnBase::name(), 
			    ScalarColumn<vector<DT> >::stringLength());
      if (nRows>0) target->v.reserve(nRows);
      for (int i=0; i<nRows; i++)
	target->v.push_back(v[i+rowStart]);
      return target;
    }
    virtual ColumnBase* copyRows(std::vector<bool>& vb, int nReserve) const {
      if (vb.size() > v.size())
	throw FTableError("too many input bools in copyRows() on column " 
			  + ColumnBase::name());
      ArrayColumn<DT>* target =
	new ArrayColumn<DT>(ColumnBase::name(), 
			    ScalarColumn<vector<DT> >::stringLength());
      if (nReserve<0) {
	nReserve = 0;
	for (std::vector<bool>::const_iterator i=vb.begin();
	     i != vb.end();
	     ++i)
	  if (*i) nReserve++;
      }
      target->v.reserve(nReserve);
      for (int i=0; i<vb.size(); i++)
	if (vb[i]) target->v.push_back(v[i]);
      return target;
    }
  };

  // Override the constructor & writing routines for Array column to check length
  // of strings in arrays of strings
  template<>
  inline ArrayColumn<string>::ArrayColumn(const vector<vector<string> >& in, string name,
				   long length_):
    Base(name, length_), dType2(FITS::FITSTypeOf<string>()),
    colCode2(FITS::ColumnCode<string>()) {
    if (stringLength()>=0)
      for (int i=0; i<in.size(); i++)
	for (int j=0; j<in[i].size(); j++)
	  checkLength(in[i][j]);
    v = in;
  }

  template<>
  long
  inline ArrayColumn<string>::writeCell(const vector<string>& value, long row) {
    if (stringLength()>=0) 
      for (int i=0; i<value.size(); i++)
	checkLength(value[i]);
    return Base::writeCell(value, row);
  }

  template <>
  long
  inline ArrayColumn<string>::writeCells(const vector<vector<string> >& values, long rowStart) {
    if (stringLength()>=0) {
      for (int i=0; i<values.size(); i++) 
	for (int j=0; j<values[i].size(); j++) 
	  checkLength(values[i][j]);
    }
    return Base::writeCells(values, rowStart);
  }
  
  template <class DT>
  class FixedArrayColumn: public ArrayColumn<DT> {
  private:
    typedef ArrayColumn<DT> Base;
    using ArrayColumn<DT>::v;
    const int width;
    // Function to pad input array to required width; exception if input too long
    void checkWidth(const vector<DT>& in) {
      if (in.size() > width)
	FormatAndThrow<FTableError>() 
	  << "Input of fixed-length array cell is too long for column "
	  << ColumnBase::name()
	  << " width " << width
	  << " input size " << in.size();
    }
  public:
    FixedArrayColumn(string name, long repeat, long length_=-1): 
      Base(name, length_), width(repeat) {}
    // Initialize with data assumes that all elements have size()<=repeat
    FixedArrayColumn(const vector<vector<DT> >& in, string name, long repeat, long length_=-1): 
      Base(name, length_), width(repeat) {
      v.resize(in.size());
      for (long i=0; i<in.size(); i++)
	writeCell(in[i], i);
    }

    virtual ColumnBase* duplicate() const {return new FixedArrayColumn(*this);}

    virtual long repeat() const {return width;}
  
    // Create new column that is a subset of this one:
    virtual ColumnBase* copyRows(long rowStart, long rowEnd) const {
      long nRows = rowEnd - rowStart;
      FixedArrayColumn<DT>* target = 
	new FixedArrayColumn<DT>(ColumnBase::name(), 
				 repeat(),
				 ScalarColumn<vector<DT> >::stringLength());
      if (nRows>0) target->v.reserve(nRows);
      for (int i=0; i<nRows; i++)
	target->v.push_back(v[i+rowStart]);
      return target;
    }
    virtual ColumnBase* copyRows(std::vector<bool>& vb, int nReserve) const {
      if (vb.size() > v.size())
	throw FTableError("too many input bools in copyRows() on column " + ColumnBase::name());
      FixedArrayColumn<DT>* target = 
	new FixedArrayColumn<DT>(ColumnBase::name(), 
				 repeat(),
				 ScalarColumn<vector<DT> >::stringLength());
      if (nReserve<0) {
	nReserve = 0;
	for (std::vector<bool>::const_iterator i=vb.begin();
	     i != vb.end();
	     ++i)
	  if (*i) nReserve++;
      }
      target->v.reserve(nReserve);
      for (int i=0; i<vb.size(); i++)
	if (vb[i]) target->v.push_back(v[i]);
      return target;
    }

    virtual long writeCell(const vector<DT>& value, long row) {
      checkWidth(value);
      if (v.size() < row+1) {
	// Fill out array with vectors of proper size
	insertRows(v.size(), row-v.size());
      }
      // Use normal variable-length array writer
      long out = Base::writeCell(value, row);
      // Fill out cell with null values as needed to meet fixed array length
      v[row].resize(width);
      return out;
    }

    // Override variable-array version with this one to make sure each row
    // is not too long, and is padded out
    virtual long writeCells(const vector<vector<DT> >& values, long rowStart=0) {
      // extend column to hold data
      if (rowStart + values.size() > v.size()) 
	insertRows(v.size(), rowStart + values.size() - v.size());
      for (int i=0; i<values.size(); i++)
	writeCell(values[i], i+rowStart);
      return v.size();
    }

    // override insert to put in proper-length vectors to cells
    // insertBeforeRow = -1 will add row(s) to the end.
    virtual long insertRows(long insertBeforeRow, long insertNumber) {
      if (insertBeforeRow < 0 || insertBeforeRow >= v.size())
	v.insert(v.end(), insertNumber, vector<DT>(width));
      else
	v.insert(v.begin()+insertBeforeRow, insertNumber, vector<DT>(width));
      return v.size();
    }

    virtual void resize(long newSize) {
      if (newSize > v.size())
	insertRows(v.size(), newSize - v.size());  //growing
      else if (newSize < v.size())
	v.resize(newSize); // shrinking
    }
  };

  class TableData {
  private:
    typedef std::map<string,ColumnBase*,stringstuff::LessNoCase> Index;
    Index columns;	// This is where the columns are stored.
    long rowReserve;	// New columns reserve at least this much space
    long rowCount;	// Keep this number and all column lengths in synch 
    bool isAltered;
    bool lock;
  public:
    TableData(long rowReserve_=0): rowReserve(rowReserve_), rowCount(0),
				   isAltered(false), lock(false)  {}
    ~TableData();
    TableData* duplicate() const;
    void copyFrom(const TableData& rhs);
      
    long nrows() const {return rowCount;}
    int ncols() const {return columns.size();}

    bool isChanged() const {return isAltered;}
    void clearChanged() {isAltered = false;}
    void touch() {isAltered = true;}
    bool isLocked() const {return lock;}
    void setLock() {lock = true;}

    // **** Make iterator classes to run over columns ***
    class iterator {
      friend class TableData;
    private:
      explicit iterator(Index::iterator i_): i(i_) {}
      Index::iterator i;
    public:
      iterator& operator++() {++i; return *this;}
      iterator& operator--() {--i; return *this;}
      bool operator==(const iterator& rhs) const {return i==rhs.i;}
      bool operator!=(const iterator& rhs) const {return i!=rhs.i;}
      ColumnBase* operator*() const {return i->second;}
    };
    class const_iterator {
      friend class TableData;
    private:
      explicit const_iterator(Index::const_iterator i_): i(i_) {}
      Index::const_iterator i;
    public:
      const_iterator(const iterator& rhs): i(rhs.i) {}
      const_iterator& operator++() {++i; return *this;}
      const_iterator& operator--() {--i; return *this;}
      bool operator==(const const_iterator& rhs) const {return i==rhs.i;}
      bool operator!=(const const_iterator& rhs) const {return i!=rhs.i;}
      const ColumnBase* operator*() const {return i->second;}
    };

    const_iterator begin() const {return const_iterator(columns.begin());}
    const_iterator end() const {return const_iterator(columns.end());}

    const_iterator const_begin() {return const_iterator(columns.begin());}
    const_iterator const_end()   {return const_iterator(columns.end());}

    iterator begin() {checkLock("begin()"); return iterator(columns.begin());}
    iterator end() {checkLock("end()"); return iterator(columns.end());}

    const ColumnBase* constColumn(string colname) const;

    const ColumnBase* operator[](string colname) const {
      return constColumn(colname);
    }
    ColumnBase* operator[](string colname);

    vector<string> listColumns() const;
    bool hasColumn(const string& colname) const;  // Case-insensitive

    // ***** Add / remove columns:

    // Add new column given the data.
    template<class T>
    void addColumn(const vector<T>& values, string columnName, 
		   long repeat=-1, long stringLength=-1 ) {
      checkLock("addColumn()");
      add(new ScalarColumn<T>(values, columnName, stringLength));
    }
    // Specialization to recognize vector<vector> as array column
    template<class T>
    void addColumn(const vector<vector<T> >& values, string columnName, 
		   long repeat=-1, long stringLength=-1) {
      checkLock();
      if (repeat<0) add(new ArrayColumn<T>(values, columnName, stringLength));
      else add(new FixedArrayColumn<T>(values, columnName, repeat, stringLength));
    }
        
    void erase(string columnName);
    void erase(iterator i);
    void clear();

    // ***** Row erase / insert:
    void eraseRows(long rowStart, long rowEnd=-1);
    void insertRows(long insertBeforeRow, long insertNumber);

    // Row and column selection.  For each, data from this instance
    // is written into *td.  Routines know what to do when td = this.

    // First selects a row range, and also limits columns to those
    // matching one of the regexp's.
    void filter(TableData* td, long rowStart, long rowEnd,
		const vector<string>& regexps) const;
    // Select rows corresponding to true element of vb:
    void filterRows(TableData* td, vector<bool>& vb) const;

    TableData* extract(long rowStart=0, long rowEnd=-1,
		       const vector<string>& regexps
		       = vector<string>(1,".*")) const;

    template <class VT>
    bool makeVector(vector<VT>& v, const string& columnName) const {
      return (*this)[columnName]->makeVector(v);
    }

    template <class Container>
    void filterRows(TableData* td, 
		    const string& columnName, const Container& keepers) const {
      const ColumnBase* cb = (*this)[columnName];
      typedef typename Container::key_type CT;
      vector<CT> data;
      if (! cb->makeVector(data))
	throw FTableError("In filterRows(), column " + cb->name() + "data type not"
			  " convertible to container");
      vector<bool> vb(data.size(), false);
      for (int i=0; i<data.size(); i++)
	if ( keepers.find(data[i]) != keepers.end()) vb[i] = true;
      filterRows(td, vb);
    }

    // Expression evaluation, convert to type T.  Use Header hh to evaluate scalars.
    template <class T>
    void evaluate(vector<T>& result,
		  const string& expression,
		  const Header* hh) const;


    // ******* Data access:
    // Access single element (copy created)
    template <class T>
    void readCell(T& value, string columnName, long row) const {
      rangeCheck(row,columnName);
      const ScalarColumn<T>* col = dynamic_cast<const ScalarColumn<T>*> ((*this)[columnName]);
      if (!col) throw FTableError("Type mismatch reading column " + columnName);
      col->readCell(value, row);
    }

    // Or a range: rowEnd is one-past-last row; -1 means go to end.
    template <class T>
    void readCells(vector<T>& values, string columnName, 
		   long rowStart=0, long rowEnd=-1) const {
      rangeCheck(rowStart,columnName);
      if (rowEnd>0) rangeCheck(rowEnd-1,columnName);
      const ScalarColumn<T>* col = dynamic_cast<const ScalarColumn<T>*> ((*this)[columnName]);
      if (!col) {
	const ScalarColumn<int>* col2 = dynamic_cast<const ScalarColumn<int>*> ((*this)[columnName]);
	const ScalarColumn<LONGLONG>* col3 = dynamic_cast<const ScalarColumn<LONGLONG>*> ((*this)[columnName]);
	const ScalarColumn<long>* col4 = dynamic_cast<const ScalarColumn<long>*> ((*this)[columnName]);
	/**/cerr << "col " << col << " int? " << col2 << " longlong " << col3 
		 << " long " << col4 <<  endl;
	/**/cerr << sizeof(int) << " " << sizeof(long) << " " << sizeof(LONGLONG) << endl;
	/**/cerr << typeid(long).name() << " LONGLONG: " << typeid(LONGLONG).name() << endl;
	/**/cerr << FITS::ColumnCode<T>() << " on cccolumn type " << (*this)[columnName]->elementType()
		 << " repeat " << (*this)[columnName]->repeat()
		 << endl;
	throw FTableError("Type mismatch reading column " + columnName);
      }
      col->readCells(values, rowStart, rowEnd);
    }

    // Writing single or multiple cells:
    template <class T>
    void writeCell(const T& value, string columnName, long row) {
      if (row<0) rangeCheck(row,columnName);
      checkLock("writeCell()");
      ScalarColumn<T>* col = dynamic_cast<ScalarColumn<T>*> ((*this)[columnName]);
      if (!col) throw FTableError("Type mismatch writing column " + columnName);
      int doneRows = col->writeCell(value, row);
      growRows(doneRows);
    }

    template <class T>
    void writeCells(const vector<T>& values, string columnName, long rowStart=0) {
      if (rowStart<0) rangeCheck(rowStart,columnName);
      checkLock("writeCells()");
      ScalarColumn<T>* col = dynamic_cast<ScalarColumn<T>*> ((*this)[columnName]);
      if (!col) throw FTableError("Type mismatch writing column " + columnName);
      int doneRows = col->writeCells(values, rowStart);
      growRows(doneRows);
    }

    // Specialization to recognize array-valued Cells
    template <class T>
    void writeCell(const vector<T>& value, string columnName, long row) {
      if (row<0) rangeCheck(row,columnName);
      checkLock("writeCell()");
      ArrayColumn<T>* col = dynamic_cast<ArrayColumn<T>*> ((*this)[columnName]);
      if (!col) throw FTableError("Type mismatch writing column " + columnName);
      int doneRows = col->writeCell(value, row);
      growRows(doneRows);
    }

    template <class T>
    void writeCells(const vector<vector<T> >& values, string columnName, long rowStart=0) {
      if (rowStart<0) rangeCheck(rowStart,columnName);
      checkLock("writeCells()");
      ArrayColumn<T>* col = dynamic_cast<ArrayColumn<T>*> ((*this)[columnName]);
      if (!col) throw FTableError("Type mismatch writing column " + columnName);
      touch();
      int doneRows = col->writeCells(values, rowStart);
      growRows(doneRows);
    }

  private:
    // No copy or assignment:
    TableData(const TableData& rhs);
    void operator=(const TableData& rhs);

    // call for anything that can alter data
    // Note not a const routine since anything that is forbidden by lock should
    // not be coded in a const routine.
    // Also anything the checks lock will alter data, so touch().
    void checkLock(const string& s="") {
      if (isLocked()) throw FTableLockError(s);
      touch();
    }
    // call for anything that can alter data
    void rangeCheck(long row, const string& s="") const {
#ifdef RANGE_CHECK
      if (row<0 || row >= nrows())
	FormatAndThrow<FTableRangeError>() << " accessing " << row << " of " << nrows()
					   << " " << s;
#endif
    }

    // Expand all columns to some new size:
    void growRows(long targetRows);
    // Add a column
    void add(ColumnBase* newColumn);
  };


} // namespace img
#endif // FTABLE2_H
