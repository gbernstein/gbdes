// Implementation code for Image class and related classes.
#include "Image.h"
#include <cstring>
#include <cstdio>
#include <cctype>

using namespace img;
using namespace std;

/////////////////////////////////////////////////////////////////////
// Routines for the underlying data structure ImageData
/////////////////////////////////////////////////////////////////////

template <class T>
void
ImageData<T>::clearChanged() {
  if (parent) 
    throw ImageError("ImageData::clearChanged() called for a subimage");
  isAltered = false;
  for (typename list<ImageData<T>*>::iterator i = children.begin();
       i != children.end();
       ++i)
    (*i)->isAltered = false;
}
 
// build or destroy the data array
template <class T>
void ImageData<T>::acquireArrays(Bounds<int> inBounds) {
  if (!inBounds) {
    // Do not allocate any space if the input bounds are undefined
    return;
  }

  bounds = inBounds;

  int	xsize, ysize;
  T *dataArray, *dptr, **rPtrs;
  xsize = bounds.getXMax() - bounds.getXMin() + 1;
  ysize = bounds.getYMax() - bounds.getYMin() + 1;
  
  dataArray  = new T[xsize*ysize]; 

  try {
    rPtrs = new T*[ysize];
  } catch (std::bad_alloc& b) {
    //Delete dataArray if row-pointers allocation fails
    delete [] dataArray;
    dataArray = nullptr;
    throw;  // Rethrow bad_alloc
  }
  rowPointers = rPtrs - bounds.getYMin();
  dptr = dataArray - bounds.getXMin();
  for (int i=bounds.getYMin(); i<=bounds.getYMax(); i++) {
    rowPointers[i] = dptr;
    dptr += xsize;
  }

  ownDataArray = true;
  ownRowPointers = true;
  isContiguous = true;
}

template <class T>
void
ImageData<T>::discardArrays() {
  if (ownDataArray) {
    //Free the data array(s)
    if (isContiguous) {
      delete [] const_location(bounds.getXMin(),bounds.getYMin());
    } else {
      T *dptr;
      for (int i=bounds.getYMin(); i<=bounds.getYMax(); i++) {
	dptr = rowPointers[i]+bounds.getXMin();
	delete [] dptr;
	rowPointers[i] = nullptr;
      }
    }
    ownDataArray = false;
  }
  if (ownRowPointers) {
    // Free the row pointer array:
    T **rptr;
    rptr = rowPointers + bounds.getYMin();
    delete [] rptr;
    rowPointers = nullptr;
    ownRowPointers = false;
  }
}

// image with unspecified data values:
template <class T>
ImageData<T>::ImageData(const Bounds<int> inBounds): rowPointers(nullptr),
						     parent(nullptr),
						     isAltered(false),
						     lock(false),
						     ownDataArray(false),
						     ownRowPointers(false) {
  acquireArrays(inBounds);
}

// image filled with a scalar:
template <class T>
ImageData<T>::ImageData(const Bounds<int> inBounds, 
			const T initValue): rowPointers(nullptr),
					    parent(nullptr),
					    isAltered(false),
					    lock(false),
					    ownDataArray(false),
					    ownRowPointers(false) {
  acquireArrays(inBounds);
  // Initialize the data
  Assert(isContiguous);
  long int i=0;
  long npix = bounds.getXMax() - bounds.getXMin() + 1;
  npix *= bounds.getYMax() - bounds.getYMin() + 1;

  for (T *initptr=location(bounds.getXMin(), bounds.getYMin()); 
       i<npix;
       i++, initptr++)
    *initptr = initValue;
}

// construct a subimage of a parent:
// Note that there is no assumption that parent is contiguous, its
// storage state could change.
template <class T>
ImageData<T>::ImageData(const Bounds<int> inBounds, 
			const ImageData<T>* _parent): 
  bounds(inBounds),
  parent(const_cast<ImageData<T>*> (_parent)),
  isAltered(_parent->isAltered),
  lock(_parent->lock),
  ownDataArray(false),
  ownRowPointers(false),
  rowPointers(_parent->rowPointers),
  isContiguous(false) {}

// Destructor:  Be sure to free all memory
template <class T>
ImageData<T>::~ImageData() {
  if (!children.empty() ) {
    // Do not throw from within destructor, just print error, quit
    cerr << "ERROR: ImageData::~ImageData for object that still has children"
	 <<endl;
    exit(1);
  }
  if (parent) parent->unlinkChild(this);
  discardArrays();
}

template <class T>
void
ImageData<T>::unlinkChild(ImageData<T>* child) const {
#ifdef _OPENMP
#pragma omp critical(imageio)
#endif
  {
    typename list<ImageData<T>*>::iterator cptr =
      find(children.begin(), children.end(), child);
    if (cptr==children.end() && !uncaught_exception()) 
      throw ImageError("ImageData::unlinkChild cannot find the child");
    children.erase(cptr);
  }
}

template <class T>
void
ImageData<T>::linkChild(ImageData<T>* child) const {
#ifdef _OPENMP
#pragma omp critical(imageio)
#endif
  {
    children.push_back(child);
  }
}

template <class T>
void
ImageData<T>::deleteSubimages() const {
#ifdef _OPENMP
#pragma omp critical(imageio)
#endif
  {
    while (!children.empty()) {
      delete *(children.begin());
      children.erase(children.begin());
    }
  }
}

// Get a subimage of this one
template <class T>
ImageData<T>* ImageData<T>::subimage(const Bounds<int> bsub) const {
  if (!bounds.includes(bsub)) 
    throw ImageError("Attempt to create subimage outside of ImageData"
		     " bounds");
  ImageData<T>* child = new ImageData<T>(bsub, this);
  linkChild(child);
  return child;
}
  
// Create a new (sub)image that is duplicate of this ones data
template <class T>
ImageData<T>* 
ImageData<T>::duplicate() const {
  ImageData<T>* dup = new ImageData<T>(bounds);
  // If the bounds are undefined and this is null image, we are done
  if (!bounds)
    return dup;
  
  Assert(dup->isContiguous);

  // Copy the data from old array to new array
  const T *inptr;
  T *dptr;
  long int xsize,ysize;
  xsize = bounds.getXMax() - bounds.getXMin() + 1;
  ysize = bounds.getYMax() - bounds.getYMin() + 1;
  
  if (isContiguous) {
    // Can be done is a single large copy if both contiguous
    inptr = const_location(bounds.getXMin(),bounds.getYMin());
    dptr = dup->location(bounds.getXMin(),bounds.getYMin());
    memmove( (void *)dptr, (void *)inptr, sizeof(T)*xsize*ysize);
  } else {
    // Do copy by rows
    for (int i=bounds.getYMin(); i<=bounds.getYMax(); i++) {
      inptr = const_location(bounds.getXMin(), i);
      dptr = dup->location(bounds.getXMin(), i);
      memmove( (void *)dptr, (void *)inptr, sizeof(T)*xsize);
    }
  }
  return dup;
}

// Change image size (flushes data), if data are under this object's
// control and there are no subimages to invalidate.
template <class T>
void
ImageData<T>::resize(const Bounds<int> newBounds) {
  if (newBounds==bounds) return;
  if (!ownDataArray || !ownRowPointers)
    throw ImageError("Cannot ImageData::resize() when data are"
		     " not owned by object");
  if (!children.empty())
    throw ImageError("Attempt to ImageData::resize() with subimages"
		     " in use.");
  if (parent)
    throw ImageError("Attempt to ImageData::resize() a subimage");
  checkLock("resize()");
  discardArrays();
  acquireArrays(newBounds);
}

// Make this ImageData be a copy of another
template <class T>
void
ImageData<T>::copyFrom(const ImageData<T>& rhs) {
  checkLock("copyFrom()");
  if (!children.empty())
    throw ImageError("Attempt to ImageData::copyFrom() with subimages"
		     " in use.");
  if (parent)
    throw ImageError("Attempt to ImageData::copyFrom() for a subimage");

  resize(rhs.getBounds());

  // Copy the data from old array to new array
  if (!bounds) return;  // but quit if there is no data to copy.
  
  T  *dptr;
  const T* inptr;
  long int xsize = bounds.getXMax() - bounds.getXMin() + 1;
  // Do copy by rows
  for (int i=bounds.getYMin(); i<=bounds.getYMax(); i++) {
    dptr = location(bounds.getXMin(), i);
    inptr = rhs.const_location(bounds.getXMin(), i);
    memmove( (void *)dptr, (void *)inptr, sizeof(T)*xsize);
  }
}

// Change origin of array, saving data
template <class T>
void
ImageData<T>::shift(int x0, int y0) {
  checkLock("shift()");
  if (!ownDataArray || !ownRowPointers)
    throw ImageError("Cannot ImageData::shift() when data are"
		     " not owned by object");
  if (!children.empty())
    throw ImageError("Attempt to ImageData::shift() with subimages"
		     " in use.");
  if (parent)
    throw ImageError("Attempt to ImageData::shift() for a subimage");

  if (!bounds)
    throw ImageError("Attempt to ImageData::shift() with undefined bounds");

  int dx = x0 - bounds.getXMin();
  int dy = y0 - bounds.getYMin();
  for (int i=bounds.getYMin(); i<=bounds.getYMax(); i++)
    rowPointers[i] -= dx;
  rowPointers -= dy;

  bounds.setXMin(bounds.getXMin() + dx);
  bounds.setXMax(bounds.getXMax() + dx);
  bounds.setYMin(bounds.getYMin() + dy);
  bounds.setYMax(bounds.getYMax() + dy);
}
  
/////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////
// Image<> class operations
/////////////////////////////////////////////////////////////////////////////

// Create subimage of given image.  
template <class T>
Image<T>
Image<T>::subimage(const Bounds<int> bsub) {
  // Make a subimage of current pixel array:
  ImageData<T>* subD=D->subimage(bsub);
  // New image gets this subimage along with the full header
  // of the parent.
  return Image(subD,H, new int(0), hcount);
}

// Create subimage of given image, this time const.  
// Notice it's necessary to create a non-const ImageData to build a new Image.
template <class T>
const Image<T>
Image<T>::subimage(const Bounds<int> bsub) const {
  // Make a subimage of current pixel array, return :
  ImageData<T>* subD=D->subimage(bsub);
  // New image gets this subimage along with the full header
  // of the parent. 
  return Image(subD,H, new int(0), hcount);
}
  
// Make a fresh copy of the image (with fresh header)
template <class T>
Image<T>
Image<T>::duplicate() const {
  ImageData<T>* dD=D->duplicate();
  Header* dH=H->duplicate();
  return Image(dD, dH);
}

// Image arithmetic operations
template <class T>
void
Image<T>::operator+=(const Image<T> rhs) {
  if (!getBounds().includes(rhs.getBounds())) 
    throw ImageError("this does not bound rhs in +=");
  transform_pixel(*this, rhs, plus<T>());
}
template <class T>
void
Image<T>::operator-=(const Image<T> rhs) {
  if (!getBounds().includes(rhs.getBounds())) 
    throw ImageError("this does not bound rhs in -=");
  transform_pixel(*this, rhs, minus<T>());
}
template <class T>
void
Image<T>::operator*=(const Image<T> rhs) {
  if (!getBounds().includes(rhs.getBounds())) 
    throw ImageError("this does not bound rhs in *=");
  transform_pixel(*this, rhs, multiplies<T>());
}
template <class T>
void
Image<T>::operator/=(const Image<T> rhs) {
  if (!getBounds().includes(rhs.getBounds())) 
    throw ImageError("this does not bound rhs in /=");
  transform_pixel(*this, rhs, divides<T>());
}

template <class T>
Image<T>
Image<T>::operator+(const Image<T> rhs) const {
  if (getBounds() != rhs.getBounds()) 
    throw ImageError("Mismatched image bound for +");
  Image<T> result(getBounds());
  transform_pixel(result, *this, rhs, plus<T>());
  return result;
}
template <class T>
Image<T>
Image<T>::operator-(const Image<T> rhs) const {
  if (getBounds() != rhs.getBounds()) 
    throw ImageError("Mismatched image bounds for -");
  Image<T> result(getBounds());
  transform_pixel(result, *this, rhs, minus<T>());
  return result;
}
template <class T>
Image<T>
Image<T>::operator*(const Image<T> rhs) const {
  if (getBounds() != rhs.getBounds()) 
    throw ImageError("Mismatched image bounds for +");
  Image<T> result(getBounds());
  transform_pixel(result, *this, rhs, multiplies<T>());
  return result;
}
template <class T>
Image<T>
Image<T>::operator/(const Image<T> rhs) const {
  if (getBounds() != rhs.getBounds()) 
    throw ImageError("Mismatched image bounds for /");
  Image<T> result(getBounds());
  transform_pixel(result, *this, rhs, divides<T>());
  return result;
}

// instantiate for expected types
template class img::Image<double>;
template class img::Image<float>;
template class img::Image<int>;
template class img::Image<short>;
template class img::ImageData<double>;
template class img::ImageData<float>;
template class img::ImageData<int>;
template class img::ImageData<short>;
