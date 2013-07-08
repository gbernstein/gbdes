//   FitsImage<T> is a class for a FITS image that is on disk.  FITS
//   files can hold multiple "Header Data Units" (HDUs), of which the
//  first (number 0 in my code, number 1 for CFITSIO) must be an image 
//  and any of the others can be an image or a FITS table.
//    FitsImage is derived from the Hdu class, which does the opening and
//  closing of FITS files and location of appropriate HDU.
//
//  When creating a FitsImage object, you give it the filename of the FITS
//  file and an optional HDU number or name.  The opening and closing of
//  the FITS file, getting proper extension, etc., will be transparent
//  to you, and any changes you make to the data will be automatically
//  written back to the disk file.  

//  FITS files are opened and closed as needed automatically.

//  On creating a FitsImage, one can specify FITS::Flags to control
//  whether the disk file is writable.  See comments in FITStypes.h for flag meanings.
//  If the image is to be read-only, 
//  an exception is thrown if it does not exist.  

//  FitsImage<T> is a template for images with data type T.  If the
//  type T does not match the native type of the FITS disk image, 
//  conversions will be done implicitly (by CFITSIO) on reading and writing.
//  Call retype() to force conversion of a (writeable) HDU to type T.

//  Once the FITSImage is created, you can access its header info via
//  Hdu::header()
//  To use image data, you create Image<T> objects via
//  FitsImage::use()  or useConst()
//   Images obtained with the use() methods are tied to the disk data;
//   any changes to these images are written back to the disk.  Destroying
//   the FitsImage before destroying all Images in use throws an exception.
//   All the Images in use() are sharing data and headers.  Destroy
//   Images when you are done using them to free up memory buffer space.
//   If you use() a read-only Image, it will be locked and attempts at
//   write access will throw exceptions.
//  FitsImage::extract() gives back an Image that is a COPY of the disk
//   data and header.  Changes to an extracted image are NOT mirrored back
//   to disk.
//  Both use() and extract() pull out the whole image by default.  Bounds
//  may be specified, you get back the intersection of our requested
//  area and the area that actually exists on disk.  Null intersection 
//  throws an exception.

//  FitsImage::flush() writes all buffers back to disk for safety.
//  FitsImage::isNull() reports whether there are any pixels in the 
//    FitsImage (by default, new extensions are created with zero dimension).
//  FitsImage::resize() changes the image dimensions, throwing an exception
//    if any images are in use() at the time or if it's not a writeable FitsImage.

// Opening the same FITS image extension with more than one
// FitsImage object may have unwanted consequences (if any writing is done)

// const-ness of  FITSImage does not have meaning; it is the read/write status
// when the FITSImage is opened that determines whether resultant Images are locked
// for writing.

// Cast your used / extracted Images to const if you want to access pixel values
// with the Image(x,y) syntax without throwing exceptions for being locked.

#ifndef FITSIMAGE_H
#define FITSIMAGE_H

//#define FITSDEBUG

#include "Hdu.h"
#include "Image.h"
#include <list>
using std::list;

namespace img {

  template <class T=float>
  class FitsImage: public FITS::Hdu {
  public:
    FitsImage(string filename, 
	      FITS::Flags f=FITS::ReadOnly,
	      int hduNumber_=1);
    FitsImage(string filename, 
	      FITS::Flags f,
	      string hduName_);
    ~FitsImage();
    
    typedef T  value_type;

    Image<T> extract(Bounds<int> b) const;
    Image<T> extract() const;	// Get full image
    Image<T> use() const;
    Image<T> use(const Bounds<int> b) const;

    // This call forces the whole image to be read into memory.
    // Will help efficiency of extracting many images instead of getting from disk.
    void load() const;
    // Dump the buffered image to save memory, if no Images are using it.
    void unload() const;
    void unload();
    
    virtual void flush() const {}   // const FitsImage will not try to write to disk
    virtual void flush() {if (isWriteable()) {flushData(); Hdu::flush();}}
    virtual bool isChanged() const {return Hdu::isChanged() || (dptr && dptr->isChanged());}

  // ???    virtual void clear();

    Bounds<int> getBounds() const {return dptr ? dptr->getBounds() : diskBounds;}
    bool isNull() const {return !getBounds();}  //no data for image?

    // Change image size - destroys data (header is kept).
    // Throws an exception if there are any extant extracted Images
    void resize(const Bounds<int> newBounds);

    // Force image data type on disk to become T.
    // Throws an exception if disk is not writeable.
    void retype();

    // Replace extension header & image with that of an entire Image
    // (If this FitsImage has an EXTNAME and the source Image does not,
    // it keeps the old extension name.)
    // Throws an exception if there are Images in use.
    void copy(const Image<T> I);

    // ??? adopt an image???

    // Convenience functions to make new FitsImage that is a COPY of
    // the specified Image (both data and header).  The FITS image is
    // NOT a mirror of the input image, it's a duplicate.
    // The ReadWrite and CreateImage flags are implicit here.
    static void writeToFITS(string filename, 
			    const Image<T> imageIn,
			    int HDUnum=1);
    static void writeToFITS(string filename, 
			    const Image<T> imageIn,
			    string HDUname);

  private:
    // hide:
    FitsImage(const FitsImage& rhs);
    void operator=(const FitsImage& rhs);

    Bounds<int> diskBounds;
    FITS::DataType nativeType;
    mutable ImageData<T>* dptr;
    mutable int *dcount;
    bool retypeOnFlush;	//Flag set if we need to change datatype to T on next flush.

    void flushData();
    void readBounds();	// Get bounds and native type from the disk file.
    // Read pixel data for requested bounds into a new ImageData
    // Exception if the bounds are not included in bounds on disk.
    ImageData<T>* readFromDisk(const Bounds<int> b) const;
    // Replace disk contents with data.  If retypeDisk=true, change
    // datatype of FITS HDU.  Disk image acquires shape of data.
    // Exception thrown if data bounds do not start at (1,1).
    void writeToDisk(const ImageData<T>* data, bool retypeDisk);
  };

}  // namespace img

#endif
