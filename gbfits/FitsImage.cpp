// FITS Image manipulation routines.
#include "FitsImage.h"

#include <typeinfo>

using namespace img;
using namespace FITS;

const int MAX_IMAGE_DIMENSIONS=10;	//maximal number of img dimensions

/////////////////////////////////////////////////////////////////
// Convenience routines to send images straight to a FITS extension
/////////////////////////////////////////////////////////////////
template <class T>
void 
FitsImage<T>::writeToFITS(const string fname, 
			  const Image<T> imageIn,
			  const int HDUnum) {
  FitsImage<T> fi(fname, ReadWrite+CreateImage, HDUnum);
  fi.copy(imageIn);
}

template <class T>
void 
FitsImage<T>::writeToFITS(const string fname, 
			  const Image<T> imageIn,
			  const string HDUname) {
  FitsImage<T> fi(fname, ReadWrite+CreateImage, HDUname);
  fi.copy(imageIn);
}

/////////////////////////////////////////////////////////////////
// Constructors/Destructors for the FitsImage objects
/////////////////////////////////////////////////////////////////
template <class T>
FitsImage<T>::FitsImage(const string filename, 
			const Flags f,
			const int hduNumber_): Hdu(filename, HDUImage, hduNumber_, f),
					       dptr(0), dcount(0),
					       retypeOnFlush(false) {
  readBounds();	
}

// Open with extension name
template <class T>
FitsImage<T>::FitsImage(const string filename, 
			const Flags f,
			const string hduName_): Hdu(filename, HDUImage, hduName_, f),
					       dptr(0), dcount(0),
					       retypeOnFlush(false) {
  readBounds();
}

template <class T>
void
FitsImage<T>::readBounds() {
  int status = moveTo();
  int naxis, bitpix;
  long naxes[MAX_IMAGE_DIMENSIONS];
  fits_get_img_param(fptr(), MAX_IMAGE_DIMENSIONS,
		     &bitpix, &naxis, naxes, &status);
  checkCFITSIO(status, "Error opening image extension in " + getFilename());
  if (naxis==2) {
    // a good 2d image
    nativeType = Bitpix_to_DataType(bitpix);
    diskBounds = Bounds<int>( int(1), naxes[0], int(1), naxes[1]);
  } else if (naxis==0) {
    // No image data, set this up as null image
    nativeType = Bitpix_to_DataType(bitpix);
    diskBounds = Bounds<int>();   //null bounds by default
  } else {
    throw FITSError("Image is not 2d: " + getFilename());
  }  
}

template <class T>
FitsImage<T>::~FitsImage() {
#ifdef FITSDEBUG
  cerr << "Destroying FitsImage " << parent.getFilename() 
       << " HDU #" << HDUnumber << endl;
#endif
  if (dptr) {
    // Make sure nothing else is using the mirrored data still:
    if (*dcount != 1)
      throwFitsOrDump("Closing FitsImage with its data still linked, file " + getFilename());
    flushData();
    delete dptr;
    delete dcount;
  }
  // Destructor for Hdu will flush header if still needed, and close files.
}

template <class T>
void
FitsImage<T>::unload() {
  if (dptr && dcount && (*dcount)==1 && !dptr->hasChildren()) {
    // Get here if there is a data structure that no one is using
    flushData();	//write any data back to disk
    delete dptr;	// Get rid of it all.
    dptr = 0;
    delete dcount;
    dcount = 0;
  }
}

// Const version will not attempt to flush data.
template <class T>
void
FitsImage<T>::unload() const {
  if (dptr && dcount && (*dcount)==1 && !dptr->hasChildren()) {
    if (dptr->isChanged())
      throw FITSError("FITSImage::unload() const called with altered data for " + getFilename());
    delete dptr;	// Get rid of it all.
    dptr = 0;
    delete dcount;
    dcount = 0;
  }
}

template <class T>
void
FitsImage<T>::load() const {
  if (dptr) return;
  dptr = readFromDisk(diskBounds);
  dcount = new int(1);
  if (!isWriteable()) dptr->setLock();
}

// Change image size - destroys data (header is kept).
template <class T>
void 
FitsImage<T>::resize(const Bounds<int> newBounds) {
  checkWriteable("resize()");
  if (dptr && dptr->hasChildren()) 
    throw FITSError("resize() for FitsImage " + getFilename() + " with subimages in use");

  if (dptr) {
    dptr->resize(newBounds);
  } else {
    // Just create image of new size to be written to disk later.
    dptr = new ImageData<T>(newBounds);
    dcount = new int(1);
  }
  if (!diskBounds) {
    // If the disk image had null size, we want to change its datatype too when we create it.
    retype();
  }
  dptr->touch();
}

// Force image on disk to be rewritten as type T.
template <class T>
void
FitsImage<T>::retype() {
  if (FITSTypeOf<T>() == nativeType) return; // Nothing to do.
  checkWriteable("retype()");
  // Read in the data if we don't have it, since we'll need to save it for retyping anyway
  load();
  retypeOnFlush = true;
  dptr->touch();
}

template <class T>
void
FitsImage<T>::copy(const Image<T> I) {
  checkWriteable("FitsImage::copy()");
  if (dptr && dptr->hasChildren())
    throw FITSError("copy() for FitsImage " + getFilename() + " with subimages in use");

  if (hptr) {
    hptr->copyFrom(*I.header());
  } else {
    hptr = I.header()->duplicate();
    hcount = new int(1);
  }
  if (dptr) {
    dptr->copyFrom(*I.data());
  } else {
    dptr = I.data()->duplicate();
    dcount = new int(1);
  }
  // Mark both as altered
  hptr->touch();
  dptr->touch();
}

template <class T>
void
FitsImage<T>::flushData() {
  if (!dptr) return;
  if (!dptr->isChanged()) return;
  checkWriteable("FitsImage::flushData()");
  writeToDisk(dptr, retypeOnFlush);
  retypeOnFlush = false;
  dptr->clearChanged();
}

// Make a new Image that is a copy of part of this image.
template <class T>
Image<T>
FitsImage<T>::extract() const {
  return extract(getBounds());
}

// Make a new Image that is a copy of part of this image.
template <class T>
Image<T>
FitsImage<T>::extract(Bounds<int> b) const {
  if (!getBounds().includes(b) && (getBounds() || b)) {
    // Problem if requested area isn't in the stored area,
    // but not a problem if both are undefined (i.e. null images)
    FormatAndThrow<FITSError>() << "Extraction bounds [" << b
				<< "] outside FitsImage bounds [" << getBounds()
				<< "]";
  }
  Header* hh = header()->duplicate();
  ImageData<T>* dd=0;
  if (!b) {
    dd = new ImageData<T>(b); // Degenerate bounds, null image
  } else if (dptr) {
    if (getBounds()==b) {
      dd = dptr->duplicate();
    } else {
      dd = dptr->subimage(b)->duplicate();
    }
  } else {
    dd = readFromDisk(b);
  }
  return Image<T>(dd,hh);
}

template <class T>
Image<T>
FitsImage<T>::use() const {
  return use(getBounds());
}

template <class T>
Image<T>
FitsImage<T>::use(const Bounds<int> b) const {
  if (!getBounds().includes(b))
    FormatAndThrow<FITSError>() << "Use bounds [" << b
				<< "] outside FitsImage bounds [" << getBounds()
				<< "]";
  Hdu::header();	// Force header loading
  load();		// Force data load
  Assert(dptr);
  // Build the table from our data and header:
  ImageData<T>* dd = 0;
  int* dc = 0;
  if (b == dptr->getBounds()) {
    // request bounds match full image - share the ImageData for use
    dd = dptr;
    dc = dcount;
  } else {
    // use a subimage of current data
    // It will inherit lock status of current data.
    dd = dptr->subimage(b);
    dc = new int(0);
  }
  return Image<T>(dd, hptr, dc, hcount);
}

template <class T>
void
FitsImage<T>::writeToDisk(const ImageData<T>* data, bool retypeDisk) {
  Assert(isWriteable());
  Bounds<int> b = data->getBounds();
  if (!b) {
    // Degenerate image, write a zero-dimensional image extension
    // of selected type
    nativeType = FITSTypeOf<T>();
    // Resize the extension
    diskBounds = b;
    int naxis(0);
    long naxes[MAX_IMAGE_DIMENSIONS];
    int bitpix = DataType_to_Bitpix(nativeType);
    int status = moveTo();
    fits_resize_img(fptr(), bitpix, naxis, naxes, &status);
    checkCFITSIO(status, "Error making 0-dim image extension in " + getFilename());
    return;
  }
    
  if (b != diskBounds || retypeDisk) {
    // Assign native type to disk file if it was empty or if commanded:
    if (retypeDisk || !diskBounds) nativeType = FITSTypeOf<T>();
    // Resize the extension
    diskBounds = b;
    int naxis(2);
    long naxes[MAX_IMAGE_DIMENSIONS];
    int bitpix = DataType_to_Bitpix(nativeType);
    if (diskBounds.getXMin() != 1
	|| diskBounds.getYMin() != 1)
      FormatAndThrow<FITSError>() << "Origin of Image writing to FITS must be (1,1), have ["
				  << b << "] for " << getFilename();
    naxes[0] = diskBounds.getXMax();
    naxes[1] = diskBounds.getYMax();
    int status = moveTo();
    fits_resize_img(fptr(), bitpix, naxis, naxes, &status);
    checkCFITSIO(status, "Error resizing image extension in " + getFilename());
  }

  // check that dimensions of ImageData match what CFITSIO expects for this extension
  int status = moveTo();
  int naxis(2);
  long naxes[MAX_IMAGE_DIMENSIONS];
  int bitpix;
  fits_get_img_param(fptr(), MAX_IMAGE_DIMENSIONS,
		     &bitpix, &naxis, naxes, &status);
  Assert( b.getXMin()==1 && b.getYMin()==1 && naxis==2
	  && b.getXMax() == naxes[0] && b.getYMax() == naxes[1]);

  if (data->contiguousData()) {
    long nelements=(b.getXMax()-b.getXMin()+1);
    nelements *= b.getYMax() - b.getYMin() + 1;
    long firstpix[2]={b.getXMin(), b.getYMin()};
    // ??? Note this is assuming that long = LONGLONG
    fits_write_pix(fptr(), FITSTypeOf<T>(),
		   firstpix, nelements, 
		   const_cast<T*> (data->const_location(b.getXMin(),b.getYMin())), 
		   &status);
  } else {
    //write to CFITSIO row by row
    long nelements=(b.getXMax()-b.getXMin()+1);
    long firstpix[2]={b.getXMin(), b.getYMin()};
      
    for (int y=b.getYMin(); y<=b.getYMax(); y++) {
      // write row to disk
      firstpix[1] = y;
      fits_write_pix(fptr(), FITSTypeOf<T>(),
		     firstpix, nelements, 
		     const_cast<T*> (data->const_location(b.getXMin(),y)), 
		     &status);
    }
  }
  checkCFITSIO(status, "Error writing image data to " + getFilename());

}

template <class T>
ImageData<T>*
FitsImage<T>::readFromDisk(const Bounds<int> b) const {
  if (!diskBounds.includes(b))
    FormatAndThrow<FITSError>() << "readFromDisk() bounds [" << b
				<< "] outside FITS image bounds [" << diskBounds
				<< "]";
  ImageData<T>* d = new ImageData<T>(b);
  Assert(d->contiguousData());
  int status = moveTo();
  long firstpix[2];
  long lastpix[2];
  long incpix[2];
  firstpix[0] = b.getXMin();
  firstpix[1] = b.getYMin();
  lastpix[0] = b.getXMax();
  lastpix[1] = b.getYMax();
  incpix[0] = 1;
  incpix[1] = 1;
  int anynul;

  fits_read_subset(fptr(), FITSTypeOf<T>(),
		   firstpix, lastpix, incpix,
		   0,
		   d->location(b.getXMin(), b.getYMin()),
		   &anynul, &status);
  checkCFITSIO(status, "Error reading image data from " + getFilename());
  d->clearChanged();
  return d;
}

template class img::FitsImage<double>;
template class img::FitsImage<float>;
template class img::FitsImage<int>;
template class img::FitsImage<short>;
