// Execute friends-of-friends object matching in 2d for catalogs, where we expect small clusters
// and fairly uniform distribution of points in a domain so that even-split cells are efficient,
// and the difference between bounding rectangle of a matched set and the true union of circles 
// is small.

#ifndef FOF_N
#define FOF_N

#include <list>
#include <set>
#include <vector>
#include <algorithm>

#include <iostream>

namespace fof {

  using std::vector;
  using std::set;
  using std::list;

  // Forward references
  template <class P, int DIM=2>
  class Cell;
  template <class P, int DIM=2>
  class Catalog;

  template <class P, int DIM=2>
  class Match: public list<const P*>  {
    friend class Cell<P,DIM>;
    friend class Catalog<P,DIM>;
  public:
    Match(const P& point, set<Cell<P,DIM>*> cellsContainingPoint):
      list<const P*>(1,&point),
      lower(point.getX()),
      upper(point.getX()),
      cells(cellsContainingPoint) {}
    bool touches(const P& point, double rad) const {
      vector<double> xp = point.getX();
      // First check whether point's rad-square would
      // touch Match's bounding rectangle:
      bool possible = true;
      for (int i=0; i<DIM && possible; i++) 
	if ( xp[i] < lower[i]-rad || xp[i] > upper[i]+rad) possible=false;
      if (!possible) return false;
      // Now check whether point is within rad of any member point:
      double radsq = rad*rad;
      for (typename list<const P*>::const_iterator i=list<const P*>::begin();
	   i != list<const P*>::end();
	   ++i) {
	double dsq=0.;
	vector<double> x2 = (*i)->getX();
	for (int j=0; j<DIM; j++)
	  dsq += (xp[j]-x2[j])*(xp[j]-x2[j]);
	if ( dsq <= radsq) return true;
      }
      return false;
    }
    void absorb(Match* rhs) {
      // Take union of bounds:
      for (int i=0; i<DIM; i++) {
	lower[i] = std::min(lower[i], rhs->lower[i]);
	upper[i] = std::max(upper[i], rhs->upper[i]);
      }
      // Suck points out of rhs:
      splice(list<const P*>::end(), *rhs);
      // Remove rhs and add this to Match list of all cells containing rhs
      for (typename set<Cell<P,DIM>*>::iterator i=rhs->cells.begin();
	   i != rhs->cells.end();
	   ++i) {
	(*i)->erase(rhs);
	(*i)->insert(this);
      }
      // And merge in rhs set of cells touched:
      cells.insert(rhs->cells.begin(), rhs->cells.end());
    }
    void add(const P& point, set<Cell<P,DIM>*> cellsContainingPoint) {
      vector<double> xp = point.getX();
      // Expand bounds:
      for (int i=0; i<DIM; i++) {
	lower[i] = std::min(lower[i], xp[i]);
	upper[i] = std::max(upper[i], xp[i]);
      }
      // add new point and its cells
      push_back(&point);
      cells.insert(cellsContainingPoint.begin(), cellsContainingPoint.end());
    }
  private:
    // Upper and lower coords of bounding rectangle of points
    vector<double> lower;
    vector<double> upper;
    // Cells that this match has points inside:
    std::set<Cell<P>*> cells;
    // Hide
    Match(const Match& rhs);
    void operator=(const Match& rhs);
  };

  template <class P, int DIM>
  class Cell {
    friend class Catalog<P,DIM>;
  public:
    Cell(vector<double>& _lower, vector<double> _upper, Catalog<P,DIM>* _own):
      lower(_lower), upper(_upper), 
      owner(_own),
      parent(0), left(0), right(0),
      splitIndex(0) {}
    // Destructor kills children so we just destroy root Cell:
    ~Cell() {
      if (left) delete left;
      if (right) delete right;
    }
    bool contains(const P& point) const {
      vector<double> xp = point.getX();
      for (int i=0; i<DIM; i++) {
	if ( xp[i] < lower[i] || xp[i] > upper[i]) return false;
      }
      return true;
    }
    void insert(Match<P,DIM>* m) {
      matches.insert(m);
      splitCheck();
    }
    void erase(Match<P,DIM>* m) {
      matches.erase(m);
    }
    void findCellsTouching(const P& point, set<Cell<P,DIM>*>& touching) {
      // Only reach this routine if touching this cell.  Add self if leaf:
      if (!left) {
	touching.insert(this);
	return;
      }
      // Otherwise descend to left and/or right child:
      double xpt = point.getX()[splitIndex];
      if (xpt >= splitValue - owner->getRadius()) right->findCellsTouching(point, touching);
      if (xpt <= splitValue + owner->getRadius()) left->findCellsTouching(point, touching);
      return;
    }
    // Split this cell if it has more than maxMatches in it:
    void splitCheck() {
      // Split cell when more than this many matches in it:
      const int MAX_MATCHES=10;
      // Don't split cells if they are already < this multiple of the match radius:
      const double MIN_SIZE_MULTIPLE=20.;
      // return if this is not a leaf, or if it's small
      if (left) return;
      if (matches.size() <= MAX_MATCHES) return;
      // Divide in 2 along axis that is not already too small
      int firstSplitIndex = splitIndex;
      while ( upper[splitIndex] - lower[splitIndex] <= MIN_SIZE_MULTIPLE*owner->getRadius()) {
	splitIndex++;
	if (splitIndex >= DIM) splitIndex=0;
	if (splitIndex==firstSplitIndex) {
	  // All axes are already at min size, no further splitting
	  return;
	}
      }
      splitValue = 0.5*(upper[splitIndex]+lower[splitIndex]);
      vector<double> upperLeft = upper;
      upperLeft[splitIndex] = splitValue;
      vector<double> lowerRight = lower;
      lowerRight[splitIndex] = splitValue;
      left = new Cell(lower, upperLeft, owner);
      right = new Cell(lowerRight, upper, owner);
      left->parent = this;
      right->parent = this;
      int childSplit = splitIndex+1;
      if (childSplit >= DIM) childSplit = 0;
      left->splitIndex = childSplit;
      right->splitIndex = childSplit;
      // Now place each Match into one or both children:
      // and tell it what new Cells it is in instead of this one
      for (typename set<Match<P,DIM>*>::iterator i=matches.begin();
      i != matches.end();
	   ++i) {
	(*i)->cells.erase(this);
	if ((*i)->lower[splitIndex] < splitValue) {
	  left->matches.insert(*i);
	  (*i)->cells.insert(left);
	}
	if ((*i)->upper[splitIndex] >= splitValue) {
	  right->matches.insert(*i);
	  (*i)->cells.insert(right);
	}
      }
      // Don't need our own match list any more:
      matches.clear();

      // Split children if needed:
      left->splitCheck();
      right->splitCheck();
  }
  private:
    Catalog<P,DIM>* owner;
    Cell<P,DIM>* parent;
    Cell<P,DIM>* left;
    Cell<P,DIM>* right;
    int splitIndex;
    double splitValue;
    vector<double> lower;
    vector<double> upper;
    set<Match<P,DIM>*> matches;
  };


  template <class P, int DIM>
  class Catalog: public set<Match<P,DIM>*> {
    public:
    typename set<Match<P,DIM>*>::iterator iterator;
    Catalog(vector<double> lower, vector<double> upper, double radius_):
      radius(radius_), root(lower,upper,this) {}
    ~Catalog() {
      // Kill all matches:
      for (typename set<Match<P,DIM>*>::iterator i= set<Match<P,DIM>*>::begin();
	   i!=set<Match<P,DIM>*>::end();
	   i++) {delete *i;}
    }
    double getRadius() const {return radius;}
    void add(const P& point) {
      // Make sure point is in the root domain?
      // I think if we don't do this, they'll be treated as if they are at the edge
      // of the domain, as far as Cell assignment, which is ok.

      // Find all cells that could touch this point:
      set<Cell<P,DIM>*> touching;
      root.findCellsTouching(point, touching);


      // Find subset of these Cells that actually contain the point
      // and make set of all Matches in all the touching Cells
      set<Cell<P,DIM>*> containing;
      set<Match<P,DIM>*> candidates;
      for (typename set<Cell<P,DIM>*>::iterator i=touching.begin();
	   i != touching.end();
	   ++i) {
	if ((*i)->contains(point)) containing.insert(*i);
	candidates.insert( (*i)->matches.begin(), (*i)->matches.end());
      }

      // Find all Matches that touch this point:
      Match<P,DIM>* primary = 0;  // the first Match that touches this point
      for (typename set<Match<P,DIM>*>::iterator i=candidates.begin();
	   i != candidates.end();
	   ++i) {
	if ((*i)->touches(point, radius)) {
	  if (!primary) {
	    // First Match found:  add Point to this, and all Cells containing it
	    primary = *i;
	    primary->add(point, containing);
	  } else {
	    // Point touches 2 Matches, which hence match each other, so combine
	    primary->absorb(*i);
	    // Get rid of the 2nd one
	    erase(*i);
	    delete *i;
	  }
	}
      }
      // Create a new Match if Point did not touch any others
      if (!primary) {
	primary = new Match<P,DIM>(point, containing);
	// Add to Catalog's index of all Matches
	insert(primary);
	// And add to every containing Cell's list of Matches:
	for (typename set<Cell<P,DIM>*>::iterator i=containing.begin();
	     i != containing.end();
	     ++i)
	  (*i)->insert(primary);
      }
    }

  private:
    double radius; // matching radius
    Cell<P,DIM> root;	// Root of the cell tree
    vector<double> lower;
    vector<double> upper;
  };



}  // end namespace
#endif // FOF2d_H
