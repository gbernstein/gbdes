// Class representing an HDU in a FITS file.
#ifndef HDU_H
#define HDU_H

// Base class for HDU's of FITS files.  The header is owned and managed by this class,
// as is making contact with CFITSIO and the FITS files.
// The Header of this HDU will be locked if it's read-only.

// When opening an HDU, note that OverwriteFile flag implies Create.
// Any Overwrite or Create implies ReadWrite also.

#include "FITS.h"
#include "Header.h"

namespace FITS {

  class Hdu {
  public:
    Hdu(const string filename,
	HDUType type_,
	int hduNumber_,
	Flags f=ReadOnly);
    Hdu(const string filename,
	HDUType type_,
	string hduName_,
	Flags f=ReadOnly);

    virtual ~Hdu();

    // Link this (writeable!) HDU to entirely new header (keeps name, version):
    void adoptHeader(img::Header* hdr, int* counter); 

    // Is this HDU allowed to be a primary extension (HDU #0)?
    bool canBePrimary() const { return (hduType==HDUImage) || (hduType==HDUAny);}
    bool isWriteable() const {return writeable;}

    string getFilename() const {return parent.getFilename();}
    int getNumber() const {return hduNumber;}
    string getName() const {return hduName;}
    HDUType getType() const {return hduType;}
    int getVersion() const {return hduVersion;}

    void setName(string hduName_);
    void setVersion(int version);

    // Move to this extension, return status
    int moveTo() const;

    virtual void flush() const {if (isWriteable()) flushHeader();}
    // isChanged() tells whether data is out of sync with disk
    virtual bool isChanged() const {return hptr ? hptr->isChanged() : false;}
    virtual void clear() { 
      checkWriteable("clear()");
      loadHeader();
      hptr->clear(); 
      hptr->touch();
    }

    // Going to let point to header out into the wild even for readonly HDU,
    // because the header itself should be locked & won't allow writing.
    img::Header* header() {loadHeader(); return hptr;}
    const img::Header* header() const {loadHeader(); return hptr;}

  protected:
    mutable img::Header* hptr;
    mutable int* hcount;  // Link counter for header
    fitsfile* fptr() const {return parent.getFitsptr();}
    // Set this when opening extension to notify derived classes that they have
    // to fill in a new extension.
    bool createdEmptyExtension;

    // Check if writeable, throw exception if not.
    void checkWriteable(string m="") const;

  private:
    // Hide copy and assignment:
    Hdu(const Hdu& rhs);
    void operator=(const Hdu& rhs);

    FitsFile parent;
    int hduNumber;
    string hduName;
    int hduVersion;
    HDUType hduType;
    const bool writeable;

    // Internals: read / flush header from file
    void loadHeader() const; // read header from file.
    void flushHeader() const; // or flush back to FITS

    // Open extension at current hduNumber (-1 to append).
    // existingType says what's there now: HDUNull = nothing there.
    void openHDU(Flags f, HDUType existingType);

    // Add an empty HDU of type emptyType so it will occupy emptyNum slot.  
    // If we are asking for num=0 and there already is one there, it
    // will be overwritten.
    // Calling with emptyType=HDUNull will create NAXES=0 image.
    void insertEmpty(int emptyType, int emptyNum, string emptyName="");
  }; // end class Hdu

} // end namespace FITS
#endif // HDU_H
