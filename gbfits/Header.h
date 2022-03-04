// File declaring the data structure that can mirror FITS headers.

#ifndef FHEADER_H
#define FHEADER_H

#include <sstream>
#include <list>
#include <algorithm>

#include "FitsTypes.h"
#include "Std.h"

/************************  Header ****************************
 * Contains auxiliary information for Images.  This includes:
 * a list of COMMENT strings
 * a list of HISTORY strings
 * a list of other keyword-indexed records, a la FITS headers.
 *
 * Header records have a keyword string, a value, and optional comment 
 * and units strings.  Individual records have a base class
 *   HdrRecordBase
 * and the derived classes are
 *   HdrRecordNull (no value)
 *   HdrRecord<T>  (value of type T).
 * Header keywords are case-insensitive and at least for FITS are
 * limited to 8 characters.
 *
 * Usually the client will not construct Headers, but always get
 * them from Images.  The most common methods will be:
 *   append("keyword",value,"comment")
 *        ...to append a new keyword/value pair to the header.
 *   replace("keyword",value,"comment")
 *        ...replaces old keyword header, or appends if no old one
 *   getValue("keyword", value)
 *        ...returns true & fills the value if keyword is in header,
 *           returns false if keyword is not already in header.
 *   addComment("comment") 
 *        ...appends a new COMMENT record.
 *   addHistory("history")  is a new HISTORY entry.
 *
 * Less frequently the client will use a list-oriented access to all
 * the header records.
 * There is an internal pointer to the header record list.
 * It is manipulated by rewind(), atEnd(), incr().  Pointer to the
 * "current" record is current().  find() moves the pointer to next
 * record that matches a keyword.
 * append() or insert() add new records at end or at current pointer.
 * You can call either of these with a keyword,value pair, and the
 * correct type for HdrRecord<T> will be inferred from the value.
 * erase() gets rid of either a certain keyword, or the current record.
 *
 * clear() flushes all records, plus the HISTORY and COMMENT lists.
 * size() is total number of records, HISTORY, and COMMENT entries.
 *
 * A copy or assignment of a Header is a deep copy (as long as
 * all the header types T are).  Header owns all the HdrRecords
 * and will delete them upon deletion of hte Header.
 *
 * Each Header keeps an "isAltered"
 * flag so one can note whether it is unchanged since its creation or
 * last call to notAltered().
 *
 * You can't erase individual comments or history entries.
 *****************************************************************/



namespace img {

  using namespace std;

  // Exception classes:
  class HeaderError: public std::runtime_error {
  public:
    HeaderError(const string m=""): 
    std::runtime_error("img::Header Error: " + m) {}
  };
  class HeaderLockError: public HeaderError {
  public:
    HeaderLockError(const string m=""): 
      HeaderError("Write access to locked data " + m) {}
  };

  //////////////////////////////////////////////////////////////////////////
  // Auxiliary information held for all images
  //////////////////////////////////////////////////////////////////////////

  // First the classes that are individual Header records:
  string KeyFormat(const string input);

  // Base class for all header entries
  class HdrRecordBase {
  public:
    HdrRecordBase(const string _kw, const string _com="",
		  const string _un=""): keyword(KeyFormat(_kw)),
					comment(_com),
					units(_un) {}
    virtual ~HdrRecordBase() {}
    virtual HdrRecordBase* duplicate() const {
      return new HdrRecordBase(*this);
    }
    bool matchesKey(const string _k) const {
      return keyword==KeyFormat(_k); 
    }
    string getComment() const {return comment;}
    void setComment(const string _c) {comment=_c;}
    string getUnits() const {return units;}
    void setUnits(const string _u) {units=_u;}
    string getKeyword() const {return keyword;}
    void setKeyword(const string _k) {keyword=KeyFormat(_k);}
    virtual void reset(const string _kw, const string _com="",
		  const string _un="") {
      keyword=_kw; comment=_com; units=_un;
    }

    // Two functions needed for easy interface to CFITSIO:
    virtual FITS::DataType dataType() const {return FITS::Tnull;}
    // return (void *) pointer to the value, if any
    virtual void* voidPtr() {return 0;}
    virtual const void* voidPtr() const {return 0;}

    // Set value for this entry from string; return true on failure
    virtual bool   setValueString(const string _v) {return false;}
    virtual string getValueString() const {return "";}
    string writeCard() const;
				
  protected:
    string keyword;
    string comment;
    string units;
  };



  class HdrRecordNull: public HdrRecordBase {
  public:
    // no additional data over base class
    HdrRecordNull(const string _kw, const string _com="",
		  const string _un=""): HdrRecordBase(_kw,_com,_un) {}
    virtual HdrRecordNull* duplicate() const {
      return new HdrRecordNull(*this);
    }
  };

  // Header Record that holds arbitary data class:
  template <class T>
  class HdrRecord: public HdrRecordBase {
  private:
    T val;
    mutable string valString;	//string representation of value
  public:
    HdrRecord(const string _kw, 
	      const T _val,
	      const string _com="",
	      const string _un=""): val(_val), 
      HdrRecordBase(_kw, _com, _un) {}
    virtual HdrRecord* duplicate() const {
      return new HdrRecord(*this);
    }

    T& Value() {return val;}
    const T& Value() const {return val;}

    void* voidPtr() {return static_cast<void *> (&val);}
    const void* voidPtr() const {return static_cast<const void *> (&val);}

    bool setValueString(const string _v) {
      istringstream iss(_v.c_str());
      string leftover;
      return !(iss >> val) || (iss >> leftover);
    };
    string getValueString() const {
      ostringstream os;
      os  << val; 
      return os.str();
    }
    FITS::DataType dataType() const {return FITS::FITSTypeOf<T>();}
  };

  //specializations for bool
  template<>
  string
  HdrRecord<bool>::getValueString() const; 

  template<>
  bool 
  HdrRecord<bool>::setValueString(const string _v);

  //and for string - enclose in quotes
  template<>
  string 
  HdrRecord<string>::getValueString() const;

  // Also specialize double since default formatting is not good;
  // need to force printing of decimal point.
  template<>
  string 
  HdrRecord<double>::getValueString() const;

  ///////////////////////////////////////////////////////////////\
  // Now the class for Header itself:
  class Header {
  private: 
    mutable std::list<HdrRecordBase*> hlist;
    mutable std::list<HdrRecordBase*>::iterator hptr;  //current record
    bool isAltered;
    bool lock;
    std::list<string> lcomment;	//Comment and History strings
    std::list<string> lhistory;
    void checkLock(const string& s="") {
      if (isLocked()) throw HeaderLockError(s);
    }
  public:
    Header(): hlist(), hptr(hlist.begin()), isAltered(false), lock(false) {}
    Header(const Header& rhs): hlist(), hptr(hlist.begin()), 
			       isAltered(false), lock(false) {
      copyFrom(rhs);
    }
    void copyFrom(const Header& rhs) {
      checkLock("copyFrom()");
      hlist.clear(); lcomment.clear(); lhistory.clear();
      std::list<HdrRecordBase*>::const_iterator rhsptr;
      for (rhsptr=rhs.hlist.begin(); rhsptr!=rhs.hlist.end(); ++rhsptr)
	hlist.push_back( (*rhsptr)->duplicate());
      hptr = hlist.begin();
      std::list<string>::const_iterator sptr;
      for (sptr=rhs.lcomment.begin(); sptr!=rhs.lcomment.end(); ++sptr)
	lcomment.push_back(*sptr);
      for (sptr=rhs.lhistory.begin(); sptr!=rhs.lhistory.end(); ++sptr)
	lhistory.push_back(*sptr);
      touch();
    }
    Header& operator=(const Header& rhs) {
      if (this==&rhs) return *this;
      checkLock("operator=");
      copyFrom(rhs); //copyFrom checks locking and touches
      isAltered = false;
      return *this;
    }
    ~Header() {
      // Unlock for destruction:
      lock = false;
      for (hptr=hlist.begin(); hptr!=hlist.end(); ++hptr)
	delete *hptr;
    }
    Header* duplicate() const {
      return new Header(*this);
    }

    // Clear all header records, plus comments & history
    void clear() {
      checkLock("clear()");
      hlist.clear();
      hptr=hlist.begin(); 
      lcomment.clear();
      lhistory.clear();
      touch();
    }
    void reset() {clear();}

    // History/Comment accessors
    const std::list<string>& comments() const {return lcomment;}
    const std::list<string>& history() const {return lhistory;}
    void addComment(const string s) {lcomment.push_back(s); touch();}
    void addHistory(const string s) {lhistory.push_back(s); touch();}

    bool isChanged() const {return isAltered;}  //changed since creation?
    void clearChanged() {isAltered=false;}	//reset altered flag
    void touch() {isAltered=true;}
    bool isLocked() const {return lock;}
    void setLock() {lock = true;}

    // Manipulate the pointer to current header record:
    void rewind() const {hptr=hlist.begin();}
    bool atEnd() const {return hptr==hlist.end();}
    int size() const {
      return hlist.size() + lcomment.size() + lhistory.size();
    }
    HdrRecordBase* current() {checkLock("current()"); touch(); return *hptr;}
    const HdrRecordBase* constCurrent() const {return *hptr;}
    void incr() const {++hptr;}
    const HdrRecordBase* current() const {return constCurrent();}

    // Append contents of another header to this one
    // Overwrite any duplicate keywords (except HISTORY and COMMENT)
    void operator+=(const Header& rhs) {
      checkLock("operator+=()");
      if (this==&rhs) return;
      for (list<HdrRecordBase*>::const_iterator rptr=rhs.hlist.begin();
	   rptr!=rhs.hlist.end();
	   ++rptr) {
	try {erase((*rptr)->getKeyword());} catch (HeaderError &i) {}
	hlist.push_back( (*rptr)->duplicate());
      }
      lcomment.insert(lcomment.end(), 
		      rhs.lcomment.begin(), 
		      rhs.lcomment.end());
      lhistory.insert(lhistory.end(), 
		      rhs.lhistory.begin(), 
		      rhs.lhistory.end());
      touch();
    }

    // Add/remove header records, by base class ptr or keyword
    void append(HdrRecordBase* record) {checkLock("append()"); hlist.push_back(record); touch();}
    void insert(HdrRecordBase* record) {
      checkLock("insert()"); hlist.insert(hptr,record); touch();
    }
    void erase() {
      checkLock("erase()");
      delete *hptr; hptr=hlist.erase(hptr); touch();
    }
    void erase(const string kw) {
      if (find(kw)) {
	checkLock("erase()");
	delete *hptr; 
	hptr=hlist.erase(hptr); 
	touch();
      } else
	throw HeaderError("Cannot find record with keyword " + kw);
    }

    template <class T>
      void append(const string keyword, const T& value, 
		  const string comment="", const string units="") {
      checkLock("append()");
      hlist.push_back(new HdrRecord<T>(keyword, value, comment,units));
      touch();
    }
    template <class T>
      void replace(const string keyword, const T& value, 
		   const string comment="", const string units="") {
      checkLock("replace()");
      try {erase(keyword);} catch (HeaderError &i) {}
      append(keyword, value, comment, units);
    }
    void appendNull(const string keyword, 
		    const string comment="") {
      checkLock("appendNull()");
      hlist.push_back(new HdrRecordNull(keyword, comment));
      touch();
    }
    template <class T>
      void insert(const string keyword, const T& value, 
		  const string comment="", const string units="") {
      checkLock("insert()");
      hlist.insert(hptr, new HdrRecord<T>(keyword, value, comment,units));
      touch();
    }
    void insertNull(const string keyword, 
		    const string comment="") {
      checkLock("insertNull()");
      hlist.insert(hptr, new HdrRecordNull(keyword, comment));
      touch();
    }
    HdrRecordBase* find(const string keyword) {
      list<HdrRecordBase*>::iterator start(hptr);
      checkLock("find()");
      touch();	// ?? note header is marked as altered just for returning
      // a non-const pointer to header record.
      for ( ; hptr!=hlist.end(); ++hptr)
	if ((*hptr)->matchesKey(keyword)) return *hptr;
      // search from beginning to starting point
      for (hptr=hlist.begin(); hptr!=hlist.end() && hptr!=start; ++hptr)
	if ((*hptr)->matchesKey(keyword)) return *hptr;
      return 0;	//null pointer if nothing found
    }
    const HdrRecordBase* findConst(const string keyword) const {
      std::list<HdrRecordBase*>::iterator start(hptr);
      for ( ; hptr!=hlist.end(); ++hptr)
	if ((*hptr)->matchesKey(keyword)) return *hptr;
      // search from beginning to starting point
      for (hptr=hlist.begin(); hptr!=hlist.end() && hptr!=start; ++hptr)
	if ((*hptr)->matchesKey(keyword)) return *hptr;
      return 0;	//null pointer if nothing found
    }
    const HdrRecordBase* find(const string keyword) const {
      return findConst(keyword);
    }
    // Get/set the value of an existing record.  Bool returns false if
    // keyword doesn't exist or does not match type of argument.
    template <class T> 
    bool getValue(const string keyword, T& outVal) const;
    template <class T> 
    bool setValue(const string keyword, const T& inVal);

  };  

  // Get header from a character stream or vice-versa:
  std::istream& operator>>(std::istream& is, Header& h);
  std::ostream& operator<<(std::ostream& os, const Header& h);

  template <class T> 
  bool 
  Header::getValue(const string keyword, T& outVal) const {
    const HdrRecordBase* b=findConst(keyword);
    if (!b) return false;
    const HdrRecord<T> *dhdr;
    dhdr = dynamic_cast<const HdrRecord<T>*> (b);
    if (!dhdr) return false;
    outVal = dhdr->Value();
    return true;
  }

  template <class T> 
  bool 
  Header::setValue(const string keyword, const T& inVal) {
    checkLock("setValue()");
    HdrRecordBase* b=find(keyword);
    if (!b) return false;
    HdrRecord<T> *dhdr;
    dhdr = dynamic_cast<HdrRecord<T>*> (b);
    if (!dhdr) return false;
    dhdr->Value() = inVal;
    touch();
    return true;
  }

} // namespace img

#endif //FHEADER_H
