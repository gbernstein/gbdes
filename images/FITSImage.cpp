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
FitsImage<T>::FitsImage(const string fname, 
			const Flags f,
			const int HDUnumber_): Hdu(filename, HDUImage, hduNumber_, f),
					       dptr(0), dcount(0),
					       retypeOnFlush(false) {
  readBounds();	
}

// Open with extension name
template <class T>
FitsImage<T>::FitsImage(const string fname, 
			const Flags f,
			const string HDUname): Hdu(filename, HDUImage, hduName_, f),
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
FitsImage<T>::unload() const {
  if (dptr && dcount && (*dcount)==1 && !dptr->hasChildren()) {
    // Get here if there is a data structure that no one is using
    flushData();	//write any data back to disk
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
  Header* hh;
  int* hc;
  ImageData<T>* dd;
  int* dc;
  I.showGuts(hh, hc, dd, dc);
  if (dptr) {
    dptr->copyFrom(*dd);
  } else {
    dptr = dd->duplicate();
    dcount = new int(1);
  }
  // Mark both as altered
  hptr->touch();
  dptr->touch();
}

template <class T>
void
FitsImage<T>::flushData() const {
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
  if (!getBounds().includes(b))
    FormatAndThrow<FITSError>() << "Extraction bounds [" << b
				<< "] outside FitsImage bounds [" << getBounds()
				<< "]";
  Header* hh = header()->duplicate();
  ImageData<T>* dd=0;
  if (dptr) {
    if (getBounds()==b) {
      dd = dptr->duplicate();
    } else {
      dd = dptr->subimage(b)->duplicate();
    }
  } else {
    dd = readFromDisk(b);
  }
  return Image<T>(hh,dd);
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
  return new ImageData<T>(hptr, hcount, dd, dc);
}

  // Fill the ImageData structure from appropriate source:
  if (get==diskBounds && !bufferBounds.includes(get)) {
    // Read entire image directly to the extracted image, no buffer
    long nelements = (diskBounds.getYMax()-diskBounds.getYMin() + 1);
    nelements *= (diskBounds.getXMax()-diskBounds.getXMin() + 1);
    
    long firstpix[2]={diskBounds.getXMin(),diskBounds.getYMin()};
    int status(0);
    // Assuming new ImageData stores data contiguously!!! ***
    fits_read_pix(parent.getFitsptr(), FITSTypeOf<T>(),
		  firstpix, nelements, NULL,
		  idata->location(get.getXMin(), get.getYMin()),
		  NULL, &status);
    if (status) throw_CFITSIO("extract() contigous read_pix on " 
			    + parent.getFilename());
  } else {
    // Build image from buffered data, copy row by row
    touch();
    bufferMustSpan(get);
    long nbytes=sizeof(T)*(get.getXMax()-get.getXMin()+1);
    for (int y=get.getYMin(); y<=get.getYMax(); ++y)
      memcpy(idata->location(get.getXMin(), y),
	     bufferLocation(get.getXMin(), y),
	     nbytes);
  }
  return Image<T>(idata,ihdr);
}

template <class T>
void
FitsImage<T>::writeToDisk(const ImageData<T>* data, bool retypeDisk) {
  Assert(isWriteable());
  if (data->getBounds() != diskBounds || retypeDisk) {
    // Resize the extension
    diskBounds = data->getBounds();
    if (retypeDisk) nativeType = FITSTypeOf<T>();
    int naxis(2);
    long naxes[MAX_IMAGE_DIMENSIONS];
    int bitpix = DataType_to_Bitpix(nativeType);
    if (data->getBounds().getXMin() != 1
	|| data->getBounds().getYMin() != 1)
      FormatAndThrow<FITSError>() << "Origin of Image writing to FITS must be (1,1), have ["
				  << data->getBounds() << "] for " << getFilename();
    naxes[0] = newBounds.getXMax();
    naxes[1] = newBounds.getYMax();
    int status = moveTo();
    fits_resize_img(fptr(), bitpix, naxis, naxes, &status);
    checkCFITSIO(status, "Error resizing image extension in " + getFilename());
  }
  // ??? Copy data onto disk
}

template <class T>
ImageType<T>*
FitsImage<T>::readFromDisk(const Bounds<int> b) const {
  if (!diskBounds.includes(b))
    FormatAndThrow<FITSError>() << "readFromDisk() bounds [" << b
				<< "] outside FITS image bounds [" << diskBounds
				<< "]";
  d = new ImageData<T>(b);
  // ??? read data from disk
  return d;
}

template class FitsImage<float>;
template class FitsImage<int>;
template class FitsImage<short>;

////////////////////////////////////////////////////////////////////////////////////////////

  if (!(parent.getFlags() & ReadWrite))
    throw FITSError("Attempt to write() to read-only FitsImage()");
  if (!(I.getBounds().includes(b) && diskBounds.includes(b)))
    throw FITSError("Attempt to write() beyond bounds of image");

  // Choose favored region to buffer.
  Bounds<int> desired=desiredBounds(b);

  long desiredSize = desired.getXMax() - desired.getXMin() + 1;
  desiredSize *= desired.getYMax() - desired.getYMin() + 1;

  int status(0);

    // If buffering this would make something too large, then think
  // about writing directly to the disk file:
  if (desiredSize > bufferTarget) {
    // Move FITS pointer to proper HDU
    fits_movabs_hdu( parent.getFitsptr(), HDUnumber, NULL, &status);
    if (!buffer) {
      // write full image directly to disk.
      // if contiguous this can be a single call:
      // isContiguous needs a check...
      if (I.data()->contiguousData() && b.getXMin()==diskBounds.getXMin()
	  && b.getXMax()==diskBounds.getXMax()) {
	long nelements=(b.getXMax()-b.getXMin()+1);
	nelements *= b.getYMax() - b.getYMin() + 1;
	long firstpix[2]={b.getXMin(), b.getYMin()};
	fits_write_pix(parent.getFitsptr(), FITSTypeOf<T>(),
		       firstpix, nelements, 
		       I.data()->location(b.getXMin(),b.getYMin()), 
		       &status);
      } else {
	//write to CFITSIO row by row
	long nelements=(b.getXMax()-b.getXMin()+1);
	long firstpix[2]={b.getXMin(), b.getYMin()};
      
	for (int y=b.getYMin(); y<=b.getYMax(); y++) {
	  // write row to disk
	  firstpix[1] = y;
	  fits_write_pix(parent.getFitsptr(), FITSTypeOf<T>(),
			 firstpix, nelements, 
			 I.data()->location(b.getXMin(),y), 
			 &status);
	}
      }
    } else {
      // write to disk those parts not in buffer
      // write the rest to present buffer
      touch();
      alterableBounds += (b & bufferBounds);	//mark as changed
      long nelements=(b.getXMax()-b.getXMin()+1);
      long nbytes=sizeof(T)*nelements;
      long firstpix[2]={b.getXMin(), b.getYMin()};
      
      for (int y=b.getYMin(); y<=b.getYMax(); y++) {
	if (bufferBounds.includes(b.getXMin(), y)) {
	  // write row to buffer
	  memcpy(bufferLocation(b.getXMin(), y),
		 I.data()->location(b.getXMin(), y),
		 nbytes);
	} else {
	  // write row to disk
	  firstpix[1] = y;
	  fits_write_pix(parent.getFitsptr(), FITSTypeOf<T>(),
			 firstpix, nelements, 
			 I.data()->location(b.getXMin(),y), 
			 &status);
	}
      }
    }
  } else {
    // Write the whole thing to a buffer
    touch();
    bufferMustSpan(b);
    alterableBounds += b;
    // note some inefficiency here, we might have just read
    // things we're going to overwrite now. ???
    // copy data to buffer row by row
    long nbytes=sizeof(T)*(b.getXMax()-b.getXMin()+1);
    for (int y=b.getYMin(); y<=b.getYMax(); ++y) {
      memcpy(bufferLocation(b.getXMin(), y),
	     I.data()->location(b.getXMin(), y),
	     nbytes);
    }
  }
  if (status) throw_CFITSIO("write() on " 
			    + parent.getFilename());
}



	fits_create_img(parent.getFitsptr(), 
			DataType_to_Bitpix(Tfloat),
			naxis, naxes,
			&status);
	if (status) throw_CFITSIO("Constructor creating image for " +
				  parent.getFilename());
  int status(0);
  int naxis(2);
  long naxes[MAX_IMAGE_DIMENSIONS];
  int bitpix = DataType_to_Bitpix(FITSTypeOf<T>());
  if (newBounds) {
    // Image is defined, make 2d
    if (newBounds.getXMin()!=1 ||
	newBounds.getYMin()!=1)
      throw FITSError("Bounds of resized FitsImage <"
		      + parent.getFilename() + 
		      "> do not start at (1,1)");
    naxes[0] = newBounds.getXMax();
    naxes[1] = newBounds.getYMax();
  } else {
    // redefine as a zero-dimensional image
    naxis=0;
  }
  // go to correct HDU
  fits_movabs_hdu(parent.getFitsptr(), HDUnumber, NULL, &status);
  parent.flush();
  fits_resize_img(parent.getFitsptr(),
		  bitpix, naxis, naxes,
		  &status);
  if (status) throw_CFITSIO("resize() on " + parent.getFilename());
  parent.flush();
  if (status) throw_CFITSIO("resize() flush on " + parent.getFilename());
   diskBounds = newBounds;
}

// Allocate space for a new buffer of size elements
template<class T>
void
FitsImage<T>::allocateBuffer(const Bounds<int> b) const {
  if (buffer) 
    throw FITSError("allocateBuffer() called when buffer!=0 for "
		    + parent.getFilename());
  if (!b) 
    throw FITSError("allocateBuffer() called with invalid bounds for "
		    + parent.getFilename());
  long size = b.getXMax() - b.getXMin() + 1;
  size *= b.getYMax() - b.getYMin() + 1;
  makeRoom(size*sizeof(T));
#ifdef FITSDEBUG
  cerr << "  About to allocate buffer for " << parent.getFilename() 
       << " HDU #" << HDUnumber
       << " with size " << size
       << endl;
#endif
  buffer = new T[size];	//catch memory failure here??? while loop?
  bufferSize = size;
  bufferRows = b.getYMax() - b.getYMin() + 1;
  alterableBounds = Bounds<int>();	//set to nil region
  totalMemoryInUse += bufferSize*sizeof(T);
  touch();
}

// Get a range of rows from disk into current data buffer.
template <class T>
void
FitsImage<T>::readRows(const int ymin, const int ymax, bool useCurrent) const{
#ifdef FITSDEBUG
  cerr << "    readRows(" << ymin << "," << ymax << ") for image "
       << parent.getFilename()
       << " HDU #" << HDUnumber 
       << endl;
#endif
  if (!buffer || (useCurrent && !bufferBounds)
      || ymin < diskBounds.getYMin()
      || ymax > diskBounds.getYMax()
      || (ymax-ymin+1) > bufferRows ) 
    throw FITSError("Bad bounds or absent buffer in readRows()");
  
  // Move the FITSFile to the proper HDU
  int status(0);
  long firstpix[2];
  long nelements;
  long xsize = diskBounds.getXMax()-diskBounds.getXMin() + 1;
  int readmin, readmax;
  firstpix[0]=diskBounds.getXMin();
  fits_movabs_hdu( parent.getFitsptr(), HDUnumber, NULL, &status);
  if (status) throw_CFITSIO("bufferEntireImage() locating HDU on " 
			    + parent.getFilename());

  // If new region will not overlap old data at all, just start over
  if (!bufferBounds 
      || ymin + bufferRows-1 < bufferBounds.getYMin()
      || ymax - bufferRows + 1 > bufferBounds.getYMax())
    useCurrent = false;

  if (!useCurrent) {
    // Ignore/Toss current data, just fill the array as desired
    firstRowOffset = 0;
    firstpix[1]=ymin;
    nelements = (ymax-ymin+1)*xsize;
    fits_read_pix(parent.getFitsptr(), FITSTypeOf<T>(),
		  firstpix, nelements, NULL,
		  buffer, NULL, &status);
    bufferBounds=Bounds<int>(diskBounds.getXMin(), diskBounds.getXMax(), 
			     ymin, ymax);
  } else {
    // Keep as much of current data in place as possible.

    // Any data to be read BELOW current buffer range?
    readmin = ymin;
    readmax = bufferBounds.getYMin()-1;
    if (readmin <= readmax) {
      Assert(readmax - readmin + 1 < bufferRows);
      int firstBuffRow = bufferBounds.getYMin() - firstRowOffset;

      if (readmin < firstBuffRow) {
	// Data to read in will go under bottom of the buffer.
	// First get the part at bottom of buffer
	firstpix[1]=firstBuffRow;
	nelements = (readmax-firstBuffRow+1)*xsize;
	if (nelements>0) fits_read_pix(parent.getFitsptr(), FITSTypeOf<T>(),
				       firstpix, nelements, NULL,
				       buffer, NULL, &status);
	// Now get the lower row range, to store at top of buffer
	firstpix[1] = readmin;
	nelements = (firstBuffRow - readmin)*xsize;
	firstRowOffset = readmin + bufferRows - firstBuffRow;
	T* target = buffer + xsize *  firstRowOffset;
	fits_read_pix(parent.getFitsptr(), FITSTypeOf<T>(),
		      firstpix, nelements, NULL,
		      target , NULL, &status);
      } else {
	// Data will fit continguously below existing data
	firstpix[1] = readmin;
	nelements = (readmax - readmin + 1)*(diskBounds.getXMax()
					     -diskBounds.getXMin() + 1);
	firstRowOffset = readmin - firstBuffRow;
	T* target = buffer + xsize*firstRowOffset;
	fits_read_pix(parent.getFitsptr(), FITSTypeOf<T>(),
		      firstpix, nelements, NULL,
		      target , NULL, &status);
      }
      // At this point, ymin is stored in firstRowOffset row of buffer.
      bufferBounds.setYMin(ymin);
      // and we may have overwritten the previous higher rows, so update
      // YMax:
      if (bufferBounds.getYMax() - bufferBounds.getYMin() > bufferRows-1)
	bufferBounds.setYMax(bufferBounds.getYMin() + bufferRows - 1);
    }
    // Any data to be read ABOVE current buffer range?
    readmin = bufferBounds.getYMax()+1;
    readmax = ymax;
    if (readmin <= readmax) {
      Assert(readmax - readmin + 1 < bufferRows);
      int lastBuffRow = bufferBounds.getYMin() - firstRowOffset
	+ bufferRows - 1;
      // put lastBuffRow in or above range to be read here
      while (lastBuffRow < readmin) lastBuffRow+=bufferRows;
      // Check for wrap around
      if (readmax > lastBuffRow && readmin <= lastBuffRow) {
	// Data to read in will wrap around the buffer.
	// First get the part at end of buffer
	firstpix[1]=readmin;
	nelements = (lastBuffRow-readmin+1)*xsize;
	T* target = buffer + xsize*
	  (bufferRows - 1 - lastBuffRow + readmin);
	fits_read_pix(parent.getFitsptr(), FITSTypeOf<T>(),
		      firstpix, nelements, NULL,
		      target, NULL, &status);
	// Now get the last row range, which will wrap around
	firstpix[1] = lastBuffRow+1;
	nelements = (readmax - lastBuffRow)*xsize;
	fits_read_pix(parent.getFitsptr(), FITSTypeOf<T>(),
		      firstpix, nelements, NULL,
		      buffer , NULL, &status);
      } else {
	// Data will fit continguously 
	firstpix[1] = readmin;
	nelements = (readmax - readmin + 1)*xsize;
	T* target;
	target = buffer + xsize*
	  (bufferRows - 1 - lastBuffRow + readmin);
	fits_read_pix(parent.getFitsptr(), FITSTypeOf<T>(),
		      firstpix, nelements, NULL,
		      target , NULL, &status);
      }
      bufferBounds.setYMax(ymax);
      // and we may have overwritten the previous lowest rows, so update
      // YMin and firstRowOffset
      if (bufferBounds.getYMax() - bufferBounds.getYMin() > bufferRows-1) {
	bufferBounds.setYMin(bufferBounds.getYMax() - bufferRows + 1);
	firstRowOffset = bufferRows - lastBuffRow 
	  + bufferBounds.getYMin() - 1;
	while (firstRowOffset < 0) firstRowOffset += bufferRows;
      }
    }
  }
  if (status) throw_CFITSIO("readRows() on " 
			    + parent.getFilename());
}

/////////////////////////////////////////////////////////////
// Write a range of rows from data buffer back to disk.
/////////////////////////////////////////////////////////////
template <class T>
void
FitsImage<T>::writeRows(const int ymin, const int ymax) const {
#ifdef FITSDEBUG
  cerr << "writeRows(" << ymin << "," << ymax << ") for image "
       << parent.getFilename() 
       << "HDU #" << HDUnumber
       << endl;
#endif
  if (!buffer 
      || ymin < diskBounds.getYMin()
      || ymax > diskBounds.getYMax()
      || ymin < bufferBounds.getYMin()
      || ymax > bufferBounds.getYMax())
    throw FITSError("Bad bounds or absent buffer in writeRows()");
  
  // Move the FITSFile to the proper HDU
  int status(0);
  long firstpix[2];
  long nelements;
  long xsize = diskBounds.getXMax() -diskBounds.getXMin() + 1;

  firstpix[0]=diskBounds.getXMin();
  fits_movabs_hdu( parent.getFitsptr(), HDUnumber, NULL, &status);
  if (status) throw_CFITSIO("writeRows() moving to HDU on " 
			    + parent.getFilename());

  // See if the region to write wraps around the buffer
  int lastBuffRow = bufferBounds.getYMin() - firstRowOffset 
    + bufferRows -1;
  if (ymin <= lastBuffRow && ymax > lastBuffRow) {
    // Need to do the writing in 2 sections.  First the lower rows:
    firstpix[1]=ymin;
    nelements = (lastBuffRow-ymin+1)*xsize;
    fits_write_pix(parent.getFitsptr(), FITSTypeOf<T>(),
		   firstpix, nelements, 
		   bufferLocation(diskBounds.getXMin(), ymin), 
		   &status);

    // Now the upper range, which wraps around to beginning of buffer
    firstpix[1] = lastBuffRow+1;
    nelements = (ymax - lastBuffRow)*xsize;
    fits_write_pix(parent.getFitsptr(), FITSTypeOf<T>(),
		   firstpix, nelements,
		   bufferLocation(diskBounds.getXMin(), lastBuffRow+1),
		   &status);
  } else {
    // Region to write is contiguous
    firstpix[1]=ymin;
    nelements = (ymax-ymin+1)*xsize;
    fits_write_pix(parent.getFitsptr(), FITSTypeOf<T>(),
		   firstpix, nelements, 
		   bufferLocation(diskBounds.getXMin(),ymin), 
		   &status);
  }
  if (status) throw_CFITSIO("writeRows() on " 
			    + parent.getFilename());
}

      T** newrpt = makeRowPointers(i->getBounds());
      i->data()->replaceRowPointers(newrpt);

// Make a new RowPointer array for some image subsection stored in buffer
template <class T>
T** 
FitsImage<T>::makeRowPointers(const Bounds<int> b) const {
  // check bounds
  if (!buffer || !b || !bufferBounds.includes(b)) {
    throw FITSError("makeRowPointers to data not in buffer");
  }
  T** rptr = new T*[b.getYMax()-b.getYMin()+1]; //??catch memory failure
  T** dptr=rptr;
  for (int y=b.getYMin(); y<=b.getYMax(); ++y, ++dptr)
    *dptr = bufferLocation(0,y);
  return rptr - b.getYMin();
}

