// Define algorithm templates which iterate over pixels in images.
// Used for example to define image arithmetic.

#ifndef IMAGE_ALGORITHM_H
#define IMAGE_ALGORITHM_H

namespace img {

  // Execute function on each pixel value
  template <class Img, class Op>
  Op for_each_pixel(Img I, Op f) {
    for (int i=I.yMin(); i<=I.yMax(); i++)
      f=for_each(I.rowBegin(i), I.rowEnd(i), f);
    return f;
  }

  // Execute function on a range of pixels
  template <class Img, class Op>
  Op for_each_pixel(Img I, Bounds<int> b, Op f) {
    if (!I.getBounds().includes(b))
      throw ImageError("for_each_pixel range exceeds image range");
    for (int i=b.getYMin(); i<=b.getYMax(); i++)
      f=for_each(I.getIterator(b.getXMin(),i), 
	       I.getIterator(b.getXMax()+1,i), f);
    return f;
  }

  // Replace image with function of itself
  template <class Img, class Op>
  Op transform_pixel(Img I, Op f) {
    for (int y=I.yMin(); y<=I.yMax(); y++) {
      typename Img::iterator ee=I.rowEnd(y);      
      for (typename Img::iterator it=I.rowBegin(y);
	   it!=ee; 
	   ++it) 
	*it=f(*it);
    }
    return f;
  }

  // Replace image with function of itself over range
  template <class Img, class Op>
  Op transform_pixel(Img I, Bounds<int> b, Op f) {
    if (!I.getBounds().includes(b))
      throw ImageError("transform_pixel range exceeds image range");
    for (int y=b.getYMin(); y<=b.getYMax(); y++) {
      typename Img::iterator ee=I.getIterator(b.getXMax()+1,y);      
      for (typename Img::iterator it=I.getIterator(b.getXMin(),y);
	   it!=ee; 
	   ++it) 
	*it=f(*it);
    }
    return f;
  }

  // Add function of pixel coords to image
  template <class Img, class Op>
  Op add_function_pixel(Img I, Op f) {
    for (int y=I.yMin(); y<=I.yMax(); y++) {
      int x=I.xMin();
      typename Img::iterator ee=I.rowEnd(y);      
      for (typename Img::iterator it=I.rowBegin(y);
	   it!=ee; 
	   ++it, ++x) 
	*it+=f(x,y);
    }
    return f;
  }

  // Add function of pixel coords to image over a range
  template <class Img, class Op>
  Op add_function_pixel(Img I, Bounds<int> b, Op f) {
    if (b && !I.getBounds().includes(b))
      throw ImageError("add_function_pixel range exceeds image range");
    for (int y=b.getYMin(); y<=b.getYMax(); y++) {
      int x=b.getXMin();
      typename Img::iterator ee=I.getIterator(b.getXMax()+1,y);      
      for (typename Img::iterator it=I.getIterator(b.getXMin(),y);
	   it!=ee; 
	   ++it, ++x) 
	*it+=f(x,y);
    }
    return f;
  }

  // Replace image with function of pixel coords
  template <class Img, class Op>
  Op fill_pixel(Img I, Op f) {
    for (int y=I.yMin(); y<=I.yMax(); y++) {
      int x=I.xMin();
      typename Img::iterator ee=I.rowEnd(y);      
      for (typename Img::iterator it=I.rowBegin(y);
	   it!=ee; 
	   ++it, ++x) 
	*it=f(x,y);
    }
    return f;
  }

  // Replace image with function of pixel coords, over specified bounds
  template <class Img, class Op>
  Op fill_pixel(Img I, Bounds<int> b, Op f) {
    if (!I.getBounds().includes(b))
      throw ImageError("add_function_pixel range exceeds image range");
    for (int y=b.getYMin(); y<=b.getYMax(); y++) {
      int x=b.getXMin();
      typename Img::iterator ee=I.getIterator(b.getXMax()+1,y);      
      for (typename Img::iterator it=I.getIterator(b.getXMin(),y);
	   it!=ee; 
	   ++it, ++x) 
	*it=f(x,y);
    }
    return f;
  }

  // Assign function of 2 images to 1st.  Pixel range comes from 2nd image
  template <class Img1, class Img2, class Op>
  Op transform_pixel(Img1 I1, const Img2 I2, Op f) {
    for (int y=I2.yMin(); y<=I2.yMax(); y++) {
      typename Img1::iterator it1=I1.getIterator(I2.xMin(),y);
      typename Img2::const_iterator ee=I2.rowEnd(y);      
      for (typename Img2::const_iterator it2=I2.rowBegin(y); 
	   it2!=ee; ++it1, ++it2) *it1=f(*it1,*it2);
    }
    return f;
  }

  // Assign function of Img2 & Img3 to Img1
  template <class Img1, class Img2, class Img3, class Op>
  Op transform_pixel(Img1 I1, const Img2 I2, const Img3 I3, Op f) {
    for (int y=I1.yMin(); y<=I1.yMax(); y++) {
      typename Img2::const_iterator it2=I2.getIterator(I1.xMin(),y);
      typename Img3::const_iterator it3=I3.getIterator(I1.xMin(),y);
      typename Img1::iterator ee=I1.rowEnd(y);      
      for (typename Img1::iterator it1=I1.rowBegin(y);
	   it1!=ee; 
	   ++it1, ++it2, ++it3) 
	*it1=f(*it2,*it3);
    }
    return f;
  }

  // Assign function of 2 images to 1st over bounds
  template <class Img1, class Img2, class Op>
  Op transform_pixel(Img1 I1, const Img2 I2, Op f, Bounds<int> b) {
    if (!I1.getBounds().includes(b) || !I2.getBounds().includes(b))
      throw ImageError("transform_pixel range exceeds image range");
    for (int y=b.getYMin(); y<=b.getYMax(); y++) {
      int x=b.getXMin();
      typename Img1::iterator ee=I1.getIterator(b.getXMax()+1,y);      
      typename Img2::const_iterator it2=I2.getIterator(b.getXMin(),y);
      for (typename Img1::iterator it1=I1.getIterator(b.getXMin(),y);
	   it1!=ee; 
	   ++it1, ++it2) 
	*it1=f(*it1,*it2);
    }
    return f;
  }



} // end namespace img
#endif  // IMAGE_ALGORITHM_H
