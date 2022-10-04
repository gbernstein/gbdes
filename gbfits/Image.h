/***********************  Image ****************************
 *  Image<T> class is a 2d array of any class T.  Data are stored 
 * in a form that allows rapid iteration along rows, and a number of
 * templates are provided that execute operations on all pixels or
 * a subset of pixels in one or more images.
 *
 * Will be linked with Image.o from Image.cpp.  Note that Image types
 * float, int, and short are instantiated at end of Image.cpp.  If you
 * need other types, you'll have to add them there.
 *
 * Image is actually a handle that contains a pointer to an Header<T>
 * structure ("header") and an ImageData<T> structure ("data").  Copy
 * and assignment semantics are that a new Image refers to same data
 * and header structures as the old one.  Link counting insures deletion
 * of the header & data structures when they become unused.  To get a
 * fresh deep copy of an Image, use Image::duplicate().
 *
 * The ImageData<T> class should never be needed by the user.  It
 * is used by FITSImage class (and maybe other disk image formats) 
 * that reads/writes Image objects from disk.  ImageData does not
 * allocate any memory when created with undefined Bounds.
 *
 * All specifications of pixel areas use Bounds<int> objects - see Bounds.h.
 * Bounds<int>(x1,x2,y1,y2) is usual constructor, Bounds<int>() creates a
 * null region.
 *
 * Image::subimage creates a new image that is contained within the original
 * image AND SHARES ITS DATA.  Deleting the parent Image before any
 * of its derived subimages throws an exception.
 *
 * Iterators Image::iter and Image::const_iter are provided to traverse rows 
 * of images.  These are only valid to traverse a row, going past end of 
 * row will give unpredictable results.  Functions rowBegin() and rowEnd() 
 * give bounds for row iteration.  getIter() gets iterator to arbitrary 
 * point.  
 * Range-checked iterators are ImageChk_iter and ImageChk_const_iter, which are 
 * also typedef'd as Image::checked_iter and checked_const_iter.  Range-checked 
 * access is via Image::at() calls.  Range-checked iterators are used for 
 * all calls if
 * #define IMAGE_BOUNDS_CHECK
; * is compiled in. 
 *
 * A const Image has read-only header and data access.
 *
 * Image constructors are:
 * Image(ncol, nrows) makes new image with origin at (1,1)
 * Image(Bounds<int>) makes new image with arbitrary row/col range
 * Image(Bounds<int>, value) makes new image w/all pixels set to value.
 *   also a constructor directly from Header & ImageData structures,
 *   which should be used only by FITSImage or routines that build Images.
 *
 * You can access image elements with (int x, int y) syntax:
 *  theImage(4,12)=15.2;
 * Many unary and binary arithmetic operations are supplied, with
 * templates provided to build others quickly.
 *
 * If an image is locked, you will get an exception if you try to acquire a
 * non-constant iterator or a reference to an element via overloaded operator(x,y).
 * To avoid this, ask for ConstIterator and access elements via value(x,y).
 * If Image instance is const, this is done automatically, whether Image is locked or not.
 ********/

#ifndef IMAGE_H
#define IMAGE_H

#include <algorithm>
#include <functional>
#include <list>
#include <sstream>
#include <typeinfo>
#include <string>

#include "Std.h"
#include "Bounds.h"
#include "Header.h"

namespace img {

  // Exception classes:
  class ImageError: public std::runtime_error {
  public: 
    ImageError(const string &m=""): std::runtime_error("Image Error: " + m) {}

  };
  class ImageBounds: public ImageError {
  public: 
    ImageBounds(const string &m=""): 
      ImageError("out of bounds " + m) {}
  };

}

#include "ImageAlgorithms.h"
#include "Image2.h"   //includes the classes not needed by end-user.

namespace img {

  template <class T=float>
  class Image {
  private:
    ImageData<T>* D;	//pixel data
    mutable int* dcount;  // link count for the data structure
    Header* H;		//the header
    mutable int* hcount;

  public:
  // Default constructor builds a null image:
    explicit Image(const Bounds<int> inBounds=Bounds<int>()): 
      D(new ImageData<T>(inBounds)),
      H(new Header()),
      dcount(new int(1)),
      hcount(new int(1)) {}
    Image(const int ncol, const int nrow):
      D(new ImageData<T>(Bounds<int>(1,ncol,1,nrow))), 
      H(new Header()),
      dcount(new int(1)),
      hcount(new int(1)) {}
    explicit Image(const Bounds<int> inBounds, const T initValue): 
      D(new ImageData<T>(inBounds, initValue)),
      H(new Header()),
      dcount(new int(1)),
      hcount(new int(1)) {}
    // Copy constructor is really a data share; lock status is shared too.
    Image(const Image &rhs): 
      D(rhs.D),
      H(rhs.H), 
      dcount(rhs.dcount),
      hcount(rhs.hcount) {(*dcount)++; (*hcount)++;}
    // Same for assignment:
    Image& operator=(const Image& rhs) {
      if (&rhs == this) return *this;
      if (D!=rhs.D) {
	if (--(*dcount)==0) {delete D; delete dcount;}
	D = rhs.D; dcount=rhs.dcount; (*dcount)++;
      }
      if (H!=rhs.H) {
	if (--(*hcount)==0) {delete H; delete hcount;}
	H = rhs.H; hcount=rhs.hcount; (*hcount)++;
      }
      return *this;
    }
    ~Image() {
      if (--(*dcount)==0) {delete D; delete dcount;}
      if (--(*hcount)==0) {delete H; delete hcount;}
    }

    // Create new image with fresh copies of data & header.  It will be unlocked, unchanged.
    Image duplicate() const;
    // New image that is subimage of this (shares pixels & header data and
    // inherits the lock / change status of this).
    Image subimage(const Bounds<int> bsub);
    const Image subimage(const Bounds<int> bsub) const ;

    ////////////////////////////////////////////////////////////
    // The following methods restructure the data array.  They will throw an exception
    // if this Image is a subimage or has extant subimages
    //
    // Make this image (or just data) be a duplicate of another's.
    // All Images that refer to same data are changed.
    void copyDataFrom(const Image& rhs) {D->copyFrom(*(rhs.D));}
    void copyFrom(const Image& rhs) {
      D->copyFrom(*(rhs.D));
      *H = *(rhs.H);
    }
    // Resize the image.
    // Note all Images sharing this ImageData will be affected.
    // Data are destroyed in the process.
    void resize(const Bounds<int> newBounds) {D->resize(newBounds);}
    //
    // Shift origin of image - same caveats apply as above
    void shift(int x0, int y0) {D->shift(x0,y0);}
    ////////////////////////////////////////////////////////////


    ////////////////////////////////////////////////////////////
    // Routines here are not meant for general usage.  They
    // are for FitsImage to be able to manipulate guts of images.

    // Create from a data and a header object: note that both will be
    // deleted when this object is deleted unless [dh]count are given.  
    Image(ImageData<T>* Din, Header* Hin,
	  int* dc=new int(0), 
	  int* hc=new int(0)): D(Din), H(Hin), 
	    dcount(dc), hcount(hc) {(*dcount)++; (*hcount)++;}

    const ImageData<T>* data() const {return D;}

    ////////////////////////////////////////////////////////////

    // Header access and some shortcuts to keywords:
    Header* header() {return H;}
    const Header* header() const {return H;}

    // Get/set the value of header records.  Bool returns false if
    // keyword doesn't exist or does not match type of argument.
    template <class U> 
    bool getHdrValue(const string keyword, U& outVal) const {
      return header()->getValue(keyword, outVal);
    }
    template <class U> 
    bool setHdrValue(const string keyword, const U& inVal) {
      return header()->setValue(keyword,inVal);
    }

    // Locking and flagging changes:
    // Note that alterations to subimages mark parents as changed.
    // Clearing change flag clears subimages.
    // locking parent locks subimages.
    // Exception for attempt to clear change flag or lock a subimage.
    bool isChanged() const {return H->isChanged() || D->isChanged();}
    void clearChanged() {D->clearChanged(); H->clearChanged();}
    bool isLocked() const {return H->isLocked() || D->isLocked();}
    // Note ***there is no unlocking of data.***  This avoids having const objects change.
    void setLock() {D->setLock(); H->setLock();}


    // Element access
#ifdef IMAGE_BOUNDS_CHECK
    // Element access is checked always
    const T operator()(int xpos, int ypos) const {
      return at(xpos,ypos);
    }
    T& operator()(int xpos, int ypos) {
      return at(xpos,ypos);
    }
    // This routine is necessary for read access to a locked, non-const Image:
    T value(int xpos, int ypos) const {
      return at(xpos, ypos);
    }
#else
    // Unchecked access
    const T operator()(const int xpos, const int ypos) const {
      return D->value(xpos,ypos);
    }
    T& operator()(const int xpos, const int ypos) {
      return (*D)(xpos,ypos);
    }
    // This routine is necessary for read access to a locked, non-const Image:
    T value(int xpos, int ypos) const {
      return D->value(xpos, ypos);
    }
#endif
    // Element access - explicitly checked
    const T at(const int xpos, const int ypos) const {
      return D->value_at(xpos,ypos);
    }
    T& at(const int xpos, const int ypos) {
      return D->at(xpos,ypos);
    }
    const T value_at(const int xpos, const int ypos) const {
      return D->value_at(xpos,ypos);
    }

    // iterators, rowBegin()/end()
    typedef Chk_iterator<T> checked_iterator;
    checked_iterator getCheckedIterator(const int x, const int y) {
      return checked_iterator(D,x,y);
    }
    typedef Chk_const_iterator<T> checked_const_iterator;
    checked_const_iterator getCheckedConstIterator(const int x, const int y) const {
      return checked_const_iterator(D,x,y);
    }
#ifdef IMAGE_BOUNDS_CHECK
    typedef checked_iterator iterator;
    typedef checked_const_iterator const_iterator;
    iterator getIterator(const int x, const int y) {
      return getCheckedIterator(x,y);}
    const_iterator getConstIterator(const int x, const int y) const {
      return getCheckedConstIterator(x,y); }
#else
    typedef T* iterator;
    typedef const T* const_iterator;
    iterator getIterator(const int x, const int y) {
      return D->location(x,y); }
    const_iterator getConstIterator(const int x, const int y) const {
      return D->const_location(x,y); }
#endif
    iterator rowBegin(int r) {return getIterator(xMin(),r);}
    const_iterator constRowBegin(int r) const {return getConstIterator(xMin(),r);}
    iterator rowEnd(int r) {return getIterator(xMax()+1,r);}
    const_iterator constRowEnd(int r) const {return getConstIterator(xMax()+1,r);}

    const_iterator getIterator(const int x, const int y) const {
      return getConstIterator(x,y); }
    const_iterator rowBegin(int r) const {return constRowBegin(r);}
    const_iterator rowEnd(int r) const {return constRowEnd(r);}

    // bounds access functions
    Bounds<int> getBounds() const {return D->getBounds();}
    int	xMin() const {return D->getBounds().getXMin();}
    int	xMax() const {return D->getBounds().getXMax();}
    int	yMin() const {return D->getBounds().getYMin();}
    int	yMax() const {return D->getBounds().getYMax();}

    // Image/scalar arithmetic operations
    void  operator+=(T x) {transform_pixel(*this, [x](T y){return y + x;});}
    void  operator-=(T x) {transform_pixel(*this, [x](T y){return y - x;});}
    void  operator*=(T x) {transform_pixel(*this, [x](T y){return y * x;});}
    void  operator/=(T x) {transform_pixel(*this, [x](T y){return y / x;});}
    void  operator-() {transform_pixel(*this, [](T y){return -y;});}

    void  operator=(const T val) {transform_pixel(*this, [val](T y){return val;});}
  
    // Image/Image arithmetic ops: rhs must be subset of this
    void  operator+=(const Image<T> rhs);
    void  operator-=(const Image<T> rhs);
    void  operator*=(const Image<T> rhs);
    void  operator/=(const Image<T> rhs);

    // Image/Image arithmetic binops: bounds must match.
    Image<T> operator+(const Image<T> rhs) const;
    Image<T> operator-(const Image<T> rhs) const;
    Image<T> operator*(const Image<T> rhs) const;
    Image<T> operator/(const Image<T> rhs) const;

    // ??? mixed-type arithmetic & conversions?
  };

} //namespace img

#endif
