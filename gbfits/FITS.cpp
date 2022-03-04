// FITSFile and FITSImage manipulations
#include "FITS.h"
#include <algorithm>
using namespace FITS;

// Static bookkeeping structures:
list<const FitsioHandle*>
FitsioHandle::openFiles;

FitsFile::HMap FitsFile::readFiles;
FitsFile::HMap FitsFile::writeFiles;

// Function to unroll the CFITSIO error stack into exception string
void
FITS::throw_CFITSIO(const string m1) {
  string m= m1 + " CFITSIO Error: ";
  char ebuff[FLEN_ERRMSG];
  while (fits_read_errmsg(ebuff)) {m+= "\n\t"; m+=ebuff;}
  // Do not throw if we are already unwinding stack from
  // another thrown exception:
  if (std::uncaught_exception()) {
    cerr << "During exception processing: " << m << endl;
  } else {
    throw FITSError(m);
  }
}

// Utility to throw exception, or if already processing exception, just print error
// Use this for things that might be thrown in destructor
// But be careful, does not interrupt control flow!!!
void
FITS::throwFitsOrDump(const string err) {
  if (std::uncaught_exception())
    cerr << "During exception processing: " << err << endl;
  else
    throw FITSError(err);
}

void
FITS::flushFitsErrors() {
  fits_clear_errmsg();
}

// FitsFile constructor:  open the file to test for existence
FitsFile::FitsFile(const string& fname, Flags f) {
  typedef std::pair<string,HCount> Entry;

  // First: be careful not to overwrite a file that is already in use,
  // either R/W or RO:
  if (f & OverwriteFile) {
    if ( readFiles.find(fname)!=readFiles.end()
	 || writeFiles.find(fname)!=writeFiles.end())
      throw FITSError("OverwriteFile specified on FITS file "
		      + fname + " that is already open.");
  }
  // Is this a duplicate filename?
  HMap* usemap = (f==FITS::ReadOnly) ? &readFiles : &writeFiles;
  HMap::iterator j = usemap->find(fname);

  if (j == usemap->end()) {
    // New file, get new handle
    ffile = new FitsioHandle(fname, f);
    usemap->insert(Entry(fname,HCount(ffile,1)));
  } else {
    // it's a duplicate, use it and increment link count:
    ffile = j->second.first;
    ++(j->second.second);
  }
}

FitsFile::~FitsFile() {
  // Find ourselves in the appropriate set of FitsHandles
  HMap* usemap = isWriteable() ? &writeFiles : &readFiles;
  HMap::iterator j = usemap->find(getFilename());
  if (j==usemap->end()) {
    ostringstream oss;
    cerr << "ERROR: Could not find file <" << getFilename()
	 << "> on FitsHandle lists in ~FitsFile()";
  }
  // Decrement link count
  --(j->second.second);
  // If links down to zero, delete the FitsHandle
  if (j->second.second <= 0) {
    delete ffile;
    usemap->erase(j);
  }
}

//////////////////////////////////////////////////////////////////
// Open/close the CFITSIO files themselves.  Automatically
// keep number of open files under some limit.
//////////////////////////////////////////////////////////////////

// FitsioHandle constructor:  open the file to test for existence
FitsioHandle::FitsioHandle(const string& fname, Flags f): 
  filename(fname), fitsptr(0) {
  int status = 0;
  writeable = !(f == FITS::ReadOnly);
  // Try to open the file - prepend "!" to filename to instruct CFITSIO to overwrite
  if (f & OverwriteFile) {
    // Stomp any existing file:
    string openName = "!" + filename;
    fits_create_file(&fitsptr, openName.c_str(), 
		     &status);
  } else {
    fits_open_file(&fitsptr, filename.c_str(), 
		   isWriteable() ? READWRITE : READONLY, 
		   &status);
  }

  if (status == FILE_NOT_OPENED) {
    if (f & FITS::Create) {
      // Presumably file is not there, try creating one
      status = 0;
      flushFitsErrors();
      fits_create_file(&fitsptr, filename.c_str(), &status);
    } else {
      // Throw special exception for failure to open
      char ebuff[256];
      fits_read_errmsg(ebuff);
      string m(ebuff);
      throw FITSCantOpen(filename, m);
    }
  }
  if (status) throw_CFITSIO("Error opening FITS file " + fname);

  // Add this freshly opened file to front of the queue:
  openFiles.push_front(this);
}

FitsioHandle::~FitsioHandle() {
  // If file is currently open:
  if (fitsptr) {
    // close file with CFITSIO
    closeFile();
    // Find and remove ourself from the list of open files
    lptr me = find(openFiles.begin(), openFiles.end(), this);
    if (me==openFiles.end()) {
      cerr << "ERROR:  Did not find open file <" + filename 
	+ "> in FitsioHandle::openFiles!" << endl;;
    } else {
      openFiles.erase(me);
    }
  }
}

// Open a file if it's not already open.  Close others if needed.
void
FitsioHandle::useThis() const {
  if (fitsptr) {
    // already open; push to front of Q for recent usage
    if (*openFiles.begin() != this) {
      lptr me = find(openFiles.begin(), openFiles.end(), this);
      openFiles.erase(me); openFiles.push_front(this);
    }
  } else {
    makeRoom();
    reopenFile();
    openFiles.push_front(this);
  }
}

void
FitsioHandle::makeRoom() const {
  static int hashReportInterval=1000;	//issue warning for file hashing.
  // set to zero to disable this reporting.
  static int hashCount=0;	//count # times need to close a file
  while (openFiles.size() >= MAX_FITS_FILES_OPEN) {
    // choose the one to close - last in queue that's open
    openFiles.back()->closeFile();
    openFiles.pop_back();
    hashCount++;
    if (hashReportInterval!=0 && (hashCount%hashReportInterval)==0)
      cerr << "WARNING: possible FitsioHandle hashing, "
	   << hashCount
	   << " files closed so far"
	   << endl;
  }
}

void
FitsioHandle::reopenFile() const {
  int status=0;
  fits_open_file(&fitsptr, filename.c_str(), 
		 isWriteable() ? READWRITE : READONLY, &status);
  if (status!=0) {
    char ebuff[256];
    fits_read_errmsg(ebuff);
    string m(ebuff);
    throw FITSCantOpen(filename, m);
  }
}

void
FitsioHandle::closeFile() const {
  if (fitsptr==0) return;
  int status=0;
#ifdef FITSDEBUG
  cerr << "fits_close_file " << filename << endl;
#endif
  fits_close_file(fitsptr, &status);
  if (status && !std::uncaught_exception()) 
  if (status) throw_CFITSIO("closeFile() on " 
			    + getFilename());
  fitsptr = 0;
}

void
FitsioHandle::flush() {
  if (fitsptr==0) return;	//No need to flush if no CFITSIO buffer
  int status=0;
  fits_flush_file(fitsptr, &status);
  if (status) throw_CFITSIO("flushing " + filename);
}

/////////////////////////////////////////////////////////////////////////
// Access FitsFile information
/////////////////////////////////////////////////////////////////////////

int
FitsFile::HDUCount() const {
  int status=0, count=0;
  fits_get_num_hdus(getFitsptr(), &count, &status);
  if (status) throw_CFITSIO("HDUCount() on " + getFilename());
  return count;
}

HDUType
FitsFile::getHDUType(const int HDUnumber) const {
  int status=0, retval;
  if (HDUnumber < 0 || HDUnumber >= HDUCount()) 
    return HDUNull;
  // CFITSIO is 1-index numbers, I am 0-indexe:
  fits_movabs_hdu(getFitsptr(), HDUnumber+1, &retval, &status);
  if (status) throw_CFITSIO("getHDUType() on " + getFilename());
  return HDUType(retval);
}

HDUType
FitsFile::getHDUType(const string HDUname, int &HDUnum) const {
  int status=0, retval;
  if (HDUCount()==0) {
    HDUnum = -1;
    return HDUNull;  // empty file
  }
  fits_movnam_hdu(getFitsptr(), ANY_HDU, 
		  const_cast<char*> (HDUname.c_str()),
		  0, &status);
  if (status==BAD_HDU_NUM) {
    //name does not exist, return HDUAny
    status = 0;
    fits_clear_errmsg();
    HDUnum = -1;
    return HDUNull;
  }
  fits_get_hdu_num(getFitsptr(), &HDUnum);
  HDUnum -= 1;	//Make zero-indexed
  fits_get_hdu_type(getFitsptr(), &retval, &status);
  if (status) throw_CFITSIO("getHDUType() by name on " + getFilename());
  return HDUType(retval);
}

HDUType
FitsFile::getHDUType(const string HDUname) const {
  int junk;
  return getHDUType(HDUname, junk);
}
