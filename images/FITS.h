// FITS.h:  manipulations of FITS files

// The FitsTypes.h file contains various enumerators and flags to
// be used with FITS files.

// The user-facing class is FitsFile.
// The class FitsioHandle represents a CFITSIO fptr and is not seen by user.
// When user creates a FitsFile, it will create a FitsioHandle or make use of
// any existing handle to the same file (although we create 2 fptr's if some
// are requesting write access and others ReadOnly).
// FitsHandle will open the file with CFITSIO, but will only keep a max number
// of CFITSIO files open at one time.  If number of FitsHandles exceeds this max,
// then only the most recently used ones are kept open.

//  The idea is that the user just opens all of them he/she might want
// and doesn't worry about it.  If there is frequent access to a large 
// number of files, then there could be hashing.  A FitsioHandle
// has fitsptr=0 if it is not currently opened.

#ifndef GBFITS_H
#define GBFITS_H

#include <list>
#include <map>
#include <cstring>
#include "Std.h"

#include "FitsTypes.h"

namespace FITS {

  // Utility function to throw FITSError and dump CFITSIO message stack.
  // Note that function will not throw while already processing an exception.
  void throw_CFITSIO(const string m="");
  // throw if status is non-zero:
  inline void checkCFITSIO(int status, const string m="") {if (status) throw_CFITSIO(m);}

  // This with throw a FITSError, but just print err to cerr if already processing exception
  void throwFitsOrDump(const string err);

  // And one to flush the error message stack
  void flushFitsErrors();
  // And clear status value:
  inline void flushFitsErrors(int& status) {status=0; flushFitsErrors();}

  // Class representing a CFITSIO fptr:
  class FitsioHandle {
  public:
    FitsioHandle(const string& fname, Flags f=ReadOnly);
    ~FitsioHandle();
    // Asking for its CFITSIO pointer reopens the file
    fitsfile* getFitsptr() const {useThis(); return fitsptr;}
    string getFilename() const {return filename;}
    bool isWriteable() const {return writeable;}
    void flush();
  private:
    const string filename;
    bool writeable;
    mutable fitsfile *fitsptr;

    // A list is kept of all open fitsio files.  Most recently
    // accessed is head of list.
    static std::list<const FitsioHandle*> openFiles;
    typedef std::list<const FitsioHandle*>::iterator lptr;

    void useThis() const;	//Make sure file is open, close others if needed
    void reopenFile() const;	//Execute the CFITSIO operations to open/close
    void closeFile() const;
    void makeRoom() const;	// Close files until we have room for a new one to open.
  };

  // User's interface:
  class FitsFile {
  public:
    FitsFile(const string& fname, Flags f=ReadOnly);
    ~FitsFile();
    fitsfile* getFitsptr() const {return ffile->getFitsptr();}
    // Implicit conversion operator to allow simple use in cfitsio calls:
    operator fitsfile*() const {return ffile->getFitsptr();}
    int HDUCount() const;
    // HDUNull is returned by the below if the HDU with chosen number or name does not exist.
    HDUType getHDUType(const int HDUnumber) const;
    HDUType getHDUType(const string HDUname) const;
    HDUType getHDUType(const string HDUname, int &HDUnum) const; //return number too
    string getFilename() const {return ffile->getFilename();}
    void flush() {ffile->flush();}
    bool isWriteable() const {return ffile->isWriteable();}
  private:
    FitsioHandle* ffile;
    int  *pcount;

    // Keep two static maps of all opened opened fitsiofile's, one for
    // writeable and one for read-only, so we can use same FITSIO handle
    // for multiple files.
    // Will keep a link count for each FitsioHandle also
    typedef std::pair<FitsioHandle*,int> HCount;
    typedef std::map<string, HCount> HMap;
    static HMap readFiles;
    static HMap writeFiles;
  };


} // end namespace FITS

#endif
