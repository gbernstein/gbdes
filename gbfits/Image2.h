// Further pieces of the Image.h header.  Here are classes that
// are not going to be used directly by the programmer.
#ifndef IMAGE2_H
#define IMAGE2_H

namespace img {
  using namespace std;

  class ImageError;
  class ImageBounds;

  ////////////////////////////////////////////////////////////////
  // The pixel data structure - never used by outside programs
  ////////////////////////////////////////////////////////////////
  template <class T=float>
  class ImageData {
  public:
    // Create:
    // image with unspecified data values:
    ImageData(const Bounds<int> inBounds) ;
    // image filled with a scalar:
    ImageData(const Bounds<int> inBounds, const T initValue) ;
    // image for which the data array has been set up by someone else:
    ImageData(const Bounds<int> inBounds, 
	      T** rptrs,
	      bool _contig=false);
  
    ~ImageData();

    // Create a new ImageData that is duplicate of this one's data
    ImageData* duplicate() const;
  
    // Manage change status - any alteration in subimage propagates to parent's change flag.
    // Clearing parent change flag clears all subimages too.
    // Cannot clear change flag of nor lock a subimage
    bool isChanged() const {return isAltered;}
    void touch() {isAltered = true; if (parent) parent->touch();}
    void clearChanged();
    bool isLocked() const {return lock;}
    void setLock() {
      if (parent) throw ImageError("Attempt to lock subimage");
      lock = true;
    }

    // Make a new ImageData that is subimage of this one.  Data will be
    // shared in memory, just pixel bounds are different.  subimages
    // a.k.a. children must be destroyed before the parent.
    ImageData<T>* subimage(const Bounds<int> bsub) const;

    // These changes to the data structure will throw exception if data array
    // is not owned or if there are subimages in use, or if this is a subimage.
    void resize(const Bounds<int> newBounds);      // Data are undefined afterwards
    // Copy in the data (and size) from another image.
    void copyFrom(const ImageData<T>& rhs);
    // move origin of coordinates, preserving data.
    void shift(int x0, int y0);

    // Element access - unchecked
    const T operator()(const int xpos, const int ypos) const {
      return *location(xpos,ypos);
    }
    T& operator()(const int xpos, const int ypos) {
      return *location(xpos,ypos);
    }
    // Explicitly read-only access:
    const T value(const int xpos, const int ypos) const {
      return *location(xpos,ypos);
    }

    // Element access - checked
    const T at(const int xpos, const int ypos) const {
      checkBounds(xpos,ypos);
      return *location(xpos,ypos);
    }
    // Explicitly read-only access:
    const T value_at(const int xpos, const int ypos) const {
      checkBounds(xpos,ypos);
      return *location(xpos,ypos);
    }
    T& at(const int xpos, const int ypos) {
      checkBounds(xpos,ypos);
      return *location(xpos,ypos);
    }

    // give pointer to a pixel in the storage array, 
    // for use by routines that buffer image data for us.
    // All write access to pixels must go through this, so this is where we check
    // for locking and update alteration flags
    T*   location(const int xpos, const int ypos) {
      checkLock("location()");
      touch();
      return *(rowPointers+ypos)+xpos;
    }
    const T* const_location(const int xpos, const int ypos) const {
      return *(rowPointers+ypos)+xpos;
    }
    const T* location(int xpos, int ypos) const {
      return const_location(xpos,ypos);
    }

    // Access functions
    Bounds<int> getBounds() const {return bounds;}
    bool contiguousData() const {return isContiguous;}

    bool hasChildren() const {return !children.empty();}

  private:
    // image which will be a subimage of a parent:
    ImageData(const Bounds<int> inBounds, 
	      const ImageData<T>* _parent);
    // Hide:
    ImageData(const ImageData &); 
    ImageData& operator=(const ImageData&);

    Bounds<int>	bounds;
    mutable T **rowPointers;	// Pointers to start of the data rows
    // Pointer is zero if this is not a subimage, else points to parent
    ImageData<T>* const parent; 
    //list of subimages of this (sub)image:
    bool isAltered;
    bool lock;
    mutable list<ImageData<T>*> children;	

    // Does this object own (i.e. have responsibility for destroying):
    bool  ownDataArray;	// the actual data array
    bool  ownRowPointers;	// the rowpointer array
    mutable bool  isContiguous;	// Set if entire image is contiguous in memory

    // class utility functions:

    // Call this for anything that is going to alter data
    void checkLock(const string& s="") {
      if (isLocked()) throw ImageError(s + " called for locked ImageData");
      touch();
    }

    void checkBounds(int x, int y) const {
      if (!bounds.includes(x,y))
	FormatAndThrow<ImageBounds>() << " at (" << x << "," << y 
				      << "), bounds " << bounds;
    }

    void acquireArrays(Bounds<int> inBounds);
    void discardArrays();
    void unlinkChild(ImageData<T>* child) const;
    void linkChild(ImageData<T>* child) const;
    // ??? Not sure if the below is needed:
    void deleteSubimages() const; // delete all living subimages

  };

  //////////////////////////////////////////////////////////////////////////
  // Checked iterator for images
  //////////////////////////////////////////////////////////////////////////
  template <class T=float>
  class Chk_iterator {
  private:
    T* ptr;
    int x;	//keep track of col number
    int xMin;
    int xMax;
  public:
    Chk_iterator(ImageData<T>* D, int x_, int y): x(x_),
						  xMin(D->getBounds().getXMin()),
						  xMax(D->getBounds().getXMax()) {
      if (y<D->getBounds().getYMin() || y>D->getBounds().getYMax())
	FormatAndThrow<ImageBounds>() << " Y= " << y
				      << " range ["  << D->getBounds().getYMin()
				      << "," << D->getBounds().getYMax() << "]";
      ptr = D->location(x,y);
    }
    T& operator*() const {
      if (x<xMin || x>xMax)
	FormatAndThrow<ImageBounds>() << " X= " << x
				      << " range [" << xMin
				      << "," << xMax << "]";
      return *ptr;
    }
    Chk_iterator operator++() {++ptr; ++x; return *this;}
    Chk_iterator operator--() {++ptr; ++x; return *this;}
    Chk_iterator operator+=(int i) {ptr+=i; x+=i; return *this;}
    Chk_iterator operator-=(int i) {ptr-=i; x-=i; return *this;}
    bool operator<(const Chk_iterator rhs) const {return ptr<rhs.ptr;}
    bool operator<=(const Chk_iterator rhs) const {return ptr<=rhs.ptr;}
    bool operator>(const Chk_iterator rhs) const {return ptr>rhs.ptr;}
    bool operator>=(const Chk_iterator rhs) const {return ptr>=rhs.ptr;}
    bool operator==(const Chk_iterator rhs) const {return ptr==rhs.ptr;}
    bool operator!=(const Chk_iterator rhs) const {return ptr!=rhs.ptr;}
  };

  template <class T=float>
  class Chk_const_iterator {
  private:
    const T* ptr;
    int x;	//keep track of col number
    int xMin;
    int xMax;
  public:
    Chk_const_iterator(const ImageData<T>* D, int x_, int y): x(x_),
							      xMin(D->getBounds().getXMin()),
							      xMax(D->getBounds().getXMax()) {
      if (y<D->getBounds().getYMin() || y>D->getBounds().getYMax())
	FormatAndThrow<ImageBounds>() << " Y= " << y
				      << " range ["  << D->getBounds().getYMin()
				      << "," << D->getBounds().getYMax() << "]";
      ptr = D->const_location(x,y);
    }
    const T& operator*() const {
      if (x<xMin || x>xMax)
	FormatAndThrow<ImageBounds>() << " X= " << x
				      << " range [" << xMin
				      << "," << xMax << "]";
      return *ptr;
    }
    Chk_const_iterator operator++() {++ptr; ++x; return *this;}
    Chk_const_iterator operator--() {++ptr; ++x; return *this;}
    Chk_const_iterator operator+=(int i) {ptr+=i; x+=i; return *this;}
    Chk_const_iterator operator-=(int i) {ptr-=i; x-=i; return *this;}
    bool operator<(const Chk_const_iterator rhs) const {return ptr<rhs.ptr;}
    bool operator<=(const Chk_const_iterator rhs) const {return ptr<=rhs.ptr;}
    bool operator>(const Chk_const_iterator rhs) const {return ptr>rhs.ptr;}
    bool operator>=(const Chk_const_iterator rhs) const {return ptr>=rhs.ptr;}
    bool operator==(const Chk_const_iterator rhs) const {return ptr==rhs.ptr;}
    bool operator!=(const Chk_const_iterator rhs) const {return ptr!=rhs.ptr;}
  };

} //namespace img
#endif //IMAGE2_H
