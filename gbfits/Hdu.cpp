// Base class for a FITS Hdu
// ??? Go through all throw's and guard against re-throws, maybe also to just die
// if multithreaded??
#include "Hdu.h"
#include <cctype>
#include <sstream>

using namespace FITS;

// Functions defined below that do the real work of header I/O to FITS:
img::Header*  readFitsHeader(fitsfile *fptr);
void  writeFitsHeader(fitsfile *fptr, img::Header *ih);
void  clearFitsHeader(fitsfile *fptr);

Hdu::Hdu(const string filename,
	 HDUType type_,
	 int hduNumber_,
	 Flags f): parent(filename, f),
		   hduNumber(hduNumber_),
		   hduName(""),
		   hduVersion(0),
		   hduType(type_),
		   writeable( !(f==ReadOnly)),
		   hptr(0),
		   hcount(0)
{
  HDUType existingType = parent.getHDUType(hduNumber);
  if (existingType != HDUNull && !(f & OverwriteHDU)) {
    // Adopt any name that this extension already has:
    int status=moveTo();
    char extname[FLEN_CARD];
    *extname=0; // empty C-string
    fits_read_key(fptr(), Tstring, "EXTNAME", extname, 0, &status);
    if (status==KEY_NO_EXIST) {
      hduName = "";
      flushFitsErrors(status);
    } else {
      hduName = extname;
    }
    checkCFITSIO(status,"Getting HDU name");
  }
  openHDU(f, existingType);
}

Hdu::Hdu(const string filename,
	 HDUType type_,
	 string hduName_,
	 Flags f): parent(filename, f),
		   hduNumber(-1),
		   hduName(hduName_),
		   hduVersion(0),
		   hduType(type_),
		   writeable( !(f==ReadOnly)),
		   hptr(0),
		   hcount(0)
{
  HDUType existingType = parent.getHDUType(hduName, hduNumber);
  openHDU(f, existingType);
}

Hdu::~Hdu() {
  if (hcount && *hcount>1) {
    // If anyone else is still using this header, there is a problem
    std::ostringstream oss;
    oss << "Header still linked when destroying HDU " << hduNumber
	<< " in file " << getFilename();
    throwFitsOrDump(oss.str());
  }
  if (isWriteable() && hptr && hptr->isChanged())
    flushHeader();
  if (hptr)   delete hptr;
  if (hcount) delete hcount;
}

void 
Hdu::adoptHeader(img::Header* hdr, int* counter) {
  checkWriteable("adoptHeader()");
  if (hcount) {
    (*hcount)--;
    if (*hcount <=0) {
      delete hptr;
      delete hcount;
    }
  }
  hptr = hdr;
  hcount = counter;
  (*hcount)++;
  hptr->touch();
}


void 
Hdu::checkWriteable(string m) const {
  if (isWriteable()) return;
  std::ostringstream oss; 
  oss << "Write action to read-only HDU #" << hduNumber
      << " name <" << hduName
      << "> in file " << getFilename()
      << "; " << m;
  throw FITSError(oss.str());
}

void
Hdu::setName(string hduName_) {
  checkWriteable("setName()");
  hduName = hduName_;
  if (hptr) hptr->touch();
}

void
Hdu::setVersion(int version) {
  checkWriteable("setVersion()");
  hduVersion = version;
  if (hptr) hptr->touch();
}

int
Hdu::moveTo() const {
  int status = 0;
  int hdutype;
  fits_movabs_hdu(fptr(), hduNumber+1, &hdutype, &status);
  return status;
}

// Connect the FITSImage disk files to Header structures
void
Hdu::loadHeader() const {
  if (hptr) return;	//already have it.

  // go to correct HDU
  int status = moveTo();
  checkCFITSIO(status, "loadHeader() for " + parent.getFilename());

  hptr = readFitsHeader(fptr());
  hcount = new int(1);
  hptr->clearChanged();	//clear the alteration flag
  // Lock header if this is read-only HDU:
  if (!isWriteable()) hptr->setLock();
}

void
Hdu::flushHeader() const {
  if (!isWriteable()) {
    // This gets called from a destructor, so may not want to throw:
    string err="flushHeader() on read-only Hdu in file " + parent.getFilename();
    throwFitsOrDump(err);
  }
  if (!hptr || !hptr->isChanged()) return;	//nothing to write back to FITS

  // Header has changed, write back to file: first clear out everything
  int status=moveTo();
  clearFitsHeader(fptr());
  // Then update the extension name and version
  if (hduName.empty()) {
    fits_delete_key(fptr(), "EXTNAME", &status);
    if (status==KEY_NO_EXIST) flushFitsErrors(status);
  } else {
    fits_update_key_str(fptr(), "EXTNAME",
			const_cast<char *> (hduName.c_str()),
			0, &status);
  }
  if (hduVersion==0)  {
    fits_delete_key(fptr(), "EXTVER", &status);
    if (status==KEY_NO_EXIST) flushFitsErrors(status);
  } else {
    fits_update_key_lng(fptr(), "EXTVER",
			(long) hduVersion, 0, &status);
  }
  checkCFITSIO(status, "flushHeader() to " + parent.getFilename());
  
  // Then append the rest of the header keywords we have:
  writeFitsHeader(fptr(), hptr);

  hptr->clearChanged();		//reset the alteration flag
}

////////////////////////////////////////////////
//  The heavy lifting for opening an HDU of a FITS file:
////////////////////////////////////////////////

void
Hdu::openHDU(Flags f, HDUType existingType) {
  createdEmptyExtension = false;
  if ( (hduNumber==0) && !canBePrimary()) {
    throw FITSError("Cannot open a FITS table at HDU number 0; file " + getFilename());
  }
  if ( existingType==HDUNull && (f & (Create + OverwriteFile))) {
    // Try to make a new extension at the end of all of them or in requested spot
    // Note that OverwriteFile will imply Create for making new extensions

    // If we've asked for appending a new extension with hduNumber=-1,
    // then target an HDU number one past end
    int nHDUs = parent.HDUCount();
    if (hduNumber < 0) hduNumber = std::max( (canBePrimary() ? 0 : 1), nHDUs);

    // Fill in 0-dimensional dummy to place our HDU in desired place.
    for (int emptyHDU=nHDUs; emptyHDU<hduNumber; emptyHDU++) {
      insertEmpty(HDUNull, emptyHDU);
    }

    // write a blank version of desired extension here:
    insertEmpty(hduType, hduNumber, hduName);

    createdEmptyExtension = true;
  } else if (existingType==HDUNull) {
    // Did not find HDU and not allowed to create it:
    FormatAndThrow<FITSError>() << "Did not find HDU number " << hduNumber
				<< " or name " << hduName
				<< " in FITS file " << getFilename();
  } else if (f & OverwriteHDU) {
    // We did find our HDU number but are instructed to overwrite it,
    // so clean it out and put in a new empty bintable
    int dummy;
    int status=0;
    fits_movabs_hdu(fptr(), hduNumber+1, &dummy, &status);
    fits_delete_hdu(fptr(), &dummy, &status);
    fits_flush_file(fptr(), &status);
    checkCFITSIO(status, "Hdu::openHDU() overwrite HDU on " + getFilename());

    // Note CFITSIO will create an empty extn #0 if you delete it.
    insertEmpty(hduType, hduNumber, hduName);
    createdEmptyExtension = true;
  } else {
    // Using pre-existing extension, no overwrite.  Make sure it is
    // of the correct type.
    if ( hduType == HDUImage && existingType != HDUImage) 
      FormatAndThrow<FITSError>() << "Existing HDU number " << hduNumber
				  << " of file " << getFilename()
				  << " is of non-image HDUType " << existingType;
    else if ( (hduType==HDUBinTable || hduType==HDUAsciiTable)
	      && !(existingType==HDUBinTable || existingType==HDUAsciiTable))
      FormatAndThrow<FITSError>() << "Existing HDU number " << hduNumber
				  << " of file " << getFilename()
				  << " is of non-table HDUType " << existingType;
    // Set type to what was found (in case we requested HDUAny)
    hduType = existingType;
    // Read extension version
    int status = moveTo();
    fits_read_key(fptr(), Tint, "EXTVER", &hduVersion, 0, &status);
    if (status==KEY_NO_EXIST) {
      hduVersion = 0;
      flushFitsErrors(status);
    } 
    checkCFITSIO(status, "Reading EXTVER");
  }
}

void
Hdu::insertEmpty(int emptyType, int emptyNum, string emptyName) {
  if (emptyNum < 0)
    FormatAndThrow<FITSError>() << "Hdu::createEmpty() with target HDU number <0: " 
				<< emptyNum << " on file " << getFilename();
  int dummy;
  int status=0;
  if (emptyType == HDUNull || emptyType == HDUImage) {
    // For case either of null or 2d image extension, will just write a 
    // 0-dimensional image here for now (image will need to set bitpix and axes 
    // in a derived class)
    if (emptyNum == 0 && parent.HDUCount() > 0) {
      // Special case because the delete_hdu call to CFITSIO will just
      // automatically replace primary extension with empty HDUImage
      fits_movabs_hdu(fptr(), 1, &dummy, &status);
      fits_delete_hdu(fptr(), &dummy, &status);
    } else {
      // insert_img a new empty image ourselves:
      // Move to previous HDU if there is one:
      if (emptyNum > 0)
	fits_movabs_hdu(fptr(), emptyNum, &dummy, &status);
      long naxes[2];
      // And insert just after:
      fits_insert_img(fptr(), BYTE_IMG, 0, naxes, &status);
    }
    // If there is a name to be added, put it in header:
    if (!emptyName.empty()) {
      fits_movabs_hdu(fptr(), emptyNum+1, &dummy, &status);
      fits_write_key_str(fptr(), "EXTNAME", const_cast<char*> (emptyName.c_str()),
			 0, &status);
    }
  } else if (emptyType==HDUBinTable) {
    // We should be creating a table now.
    if (emptyNum==0) 
      throw FITSError("Attempt to create new table at HDU #0 in file"
		      + getFilename());
    fits_movabs_hdu(fptr(), emptyNum, &dummy, &status);
    char* tDummy[1];
    tDummy[0] = 0;
    fits_create_tbl(parent, BINARY_TBL, (LONGLONG) 0, 0, 
		    tDummy, tDummy, tDummy,
		    const_cast<char*> (emptyName.c_str()),
		    &status);
  } else {
    FormatAndThrow<FITSError>() << "Creation of HDU Type "
				<< emptyType
				<< " is not implemented (e.g. ASCII table)";
  }

  // Make sure this change propagates to all other extensions in CFITSIO:
  fits_flush_file(fptr(), &status);
  // Make sure we can move to this
  fits_movabs_hdu(parent, emptyNum+1, &dummy, &status);
  checkCFITSIO(status, "Creating new extension in " + getFilename());
}


/////////////////////////////////////////////////////////////////
// Reading & Writing FITS headers to Header structures
/////////////////////////////////////////////////////////////////

// These are the FITS Header keywords that specify the extensions
// and the image/table size.  They will NOT be saved into our Header
// structure, as they are read/written automatically by CFITSIO when maintaining
// the HDUs and images.
// **Keyword only has to match the letters given here to be excluded, i.e.
// ** appending a wildcard after each of these.

bool isSpecialKeyword(const string keyword) {
  static list<string> specialKeys;
  if (specialKeys.empty()) {
    // Initial list:
    specialKeys.push_back("SIMPLE");
    specialKeys.push_back("BITPIX");
    specialKeys.push_back("NAXIS");
    specialKeys.push_back("EXTEND");
    specialKeys.push_back("XTENSION");
    specialKeys.push_back("PCOUNT");
    specialKeys.push_back("GCOUNT");
    specialKeys.push_back("TFIELDS");
    specialKeys.push_back("TTYPE");
    specialKeys.push_back("TFORM");
    specialKeys.push_back("TBCOL");
    // ** The following should be tracked in FTable column structures??
    // ** but I'm just killing them so they don't get
    // ** confused when column numbers change:
    specialKeys.push_back("TDIM");
    specialKeys.push_back("TUNIT");
    specialKeys.push_back("TSCAL");
    specialKeys.push_back("TZERO");
    specialKeys.push_back("TDISP");
    specialKeys.push_back("TNULL");
    // These will be kept track of as Hdu class members:
    specialKeys.push_back("EXTNAME");
    specialKeys.push_back("EXTVER");
  }
  for (list<string>::const_iterator i=specialKeys.begin();
       i!=specialKeys.end();
       ++i) {
    if (keyword.size() < i->size()) continue;
    int j=0;
    for ( ; j<i->size(); j++)
      if (std::toupper(keyword[j])!=std::toupper((*i)[j])) break;
    if ( j >= i->size()) 
      return true; // match here!
  }
  return false;
}

// read a full header from FITS extension (assumed current HDU)

img::Header*
readFitsHeader(fitsfile *fptr) {
  img::Header* ih=new img::Header();

  int status(0);
  int nkeys;
  fits_get_hdrspace(fptr, &nkeys, NULL, &status);

  char keyword[FLEN_CARD];
  char comment[FLEN_CARD];
  char value[FLEN_CARD];
  char units[FLEN_CARD];
  char vtype;
 
  for (int ikey=1; ikey<=nkeys; ikey++) {
    fits_read_keyn(fptr,ikey,keyword,value,comment,&status);
    if (strlen(value)>0) {
      fits_get_keytype(value, &vtype, &status);
      fits_read_key_unit(fptr,keyword,units,&status);
    } else {
      vtype='N';
      units[0]=0;
    }
    checkCFITSIO(status, "readFitsHeader collecting all keys");
    // Make the HdrRecord: 
    img::HdrRecordBase* hh;
    bool badstring(false);
    string vstring=value;
    switch (vtype) {
    case 'N': 
      hh = new img::HdrRecordNull(keyword, comment);
      break;
    case 'C':
      if (vstring[0]!='\'' || vstring[vstring.size()-1]!='\'')
	cerr << "no quotes on string: [" << vstring << "]" << endl;
      hh = new img::HdrRecord<string>(keyword,
				 vstring.substr(1,vstring.size()-2),
				 comment,
				 units);
      break;
    case 'L':
      hh = new img::HdrRecord<bool>(keyword, false, comment, units);
      badstring = hh->setValueString(vstring);
      break;
    case 'I':
      hh = new img::HdrRecord<int>(keyword, 0, comment, units);
      badstring = hh->setValueString(vstring);
      break;
    case 'F':
      hh = new img::HdrRecord<double>(keyword, 0., comment, units);
      badstring = hh->setValueString(vstring);
      break;
    case 'X':
      hh = new img::HdrRecord<complex<double> >(keyword, 
					  complex<double>(), comment, units);
      badstring = hh->setValueString(vstring);
      break;
    default:
      throw img::HeaderError("Header vstring [" + vstring +
			     "] of unknown type " + vtype);
    }
    if (badstring)
      throw img::HeaderError("Header vstring [" + vstring +
			     "] could not be interpeted as type " + vtype);
    // Skip keywords that are part of FITS extension/image definitions
    if (isSpecialKeyword(hh->getKeyword())) {
      delete hh;
    } else if (hh->matchesKey("COMMENT")) {
      ih->addComment(hh->getComment());
      delete hh;
    } else if (hh->matchesKey("HISTORY")) {
      ih->addHistory(hh->getComment());
      delete hh;
    } else {
      ih->append(hh);
    }
  }
  return ih;
}

// Append contents of Header to header of a FITS file's
// current HDU:
void
writeFitsHeader(fitsfile *fptr, img::Header *ih) {
  int status = 0;
  char keyword[FLEN_CARD];
  char comment[FLEN_CARD];
  char vstring[FLEN_CARD];
  char units[FLEN_CARD];

  // First write all the Records (stop if FITS error)
  for (ih->rewind(); !status && !ih->atEnd(); ih->incr()) {
    const img::HdrRecordBase* hh=ih->constCurrent();
    DataType t=hh->dataType();
    strncpy(keyword,hh->getKeyword().c_str(),sizeof(keyword));
    strncpy(comment,hh->getComment().c_str(),sizeof(comment));
    strncpy(units,hh->getUnits().c_str(),sizeof(units));

    switch (t) {
    case FITS::Tnull:
      fits_write_key_null(fptr, keyword, comment, &status);
      break;
    case FITS::Tstring:
      { 
	const img::HdrRecord<string>* hs 
	  = dynamic_cast<const img::HdrRecord<string>*> (hh);
	strncpy(vstring,hs->Value().c_str(),sizeof(vstring));
	vstring[sizeof(vstring)-1]=0;
      }
      fits_write_key_str(fptr, keyword, vstring, comment, &status);
      if (!hh->getUnits().empty())
	fits_write_key_unit(fptr, keyword, units, &status);
      break;
    case FITS::Tlogical:
      {
	// Specialize because CFITSIO wants pointer to int for logical value:
	const img::HdrRecord<bool>* hb 
	  = dynamic_cast<const img::HdrRecord<bool>*> (hh);
	int i = hb->Value() ? 1 : 0;
	fits_write_key(fptr, t, keyword, &i, comment, &status);
	if (!hh->getUnits().empty())
	  fits_write_key_unit(fptr, keyword, units, &status);
      }
      break;
    default:
      {
	// nasty cast needed for CFITSIO:
	void *vv = const_cast<void *> (hh->voidPtr());
	fits_write_key(fptr, t, keyword, vv, comment, &status);
      }
      if (!hh->getUnits().empty())
	fits_write_key_unit(fptr, keyword, units, &status);
    }

  }
  // Then write all the COMMENT fields
  list<string>::const_iterator sptr;
  for (sptr=ih->comments().begin(); sptr!=ih->comments().end(); ++sptr) {
    strncpy(comment, sptr->c_str(), sizeof(comment));
    comment[sizeof(comment)-1]=0;
    fits_write_comment(fptr, comment, &status);
  }
  // And HISTORY strings
  for (sptr=ih->history().begin(); sptr!=ih->history().end(); ++sptr) {
    strncpy(comment, sptr->c_str(), sizeof(comment));
    comment[sizeof(comment)-1]=0;
    fits_write_history(fptr, comment, &status);
  }
  checkCFITSIO(status, "writeFitsHeader()");
}

// Clear a FITS header of all optional keywords
void
clearFitsHeader(fitsfile *fptr) {
  int status=0;
  int nkeys, nextkey;
  char keyword[FLEN_CARD];
  char comment[FLEN_CARD];
  char value[FLEN_CARD];
  char units[FLEN_CARD];
  
  // Count keys
  fits_get_hdrpos(fptr, &nkeys, &nextkey, &status);

  for (int ikey=1; !status && ikey<=nkeys; ) {
    fits_read_keyn(fptr,ikey,keyword,value,comment,&status);
    string key(keyword);
    if (!isSpecialKeyword(key)) {
      fits_delete_record(fptr, ikey, &status);
      fits_get_hdrpos(fptr, &nkeys, &nextkey, &status);
    } else {
      ikey++;
    }
  }
  checkCFITSIO(status, "clearFitsHeader()");
}

