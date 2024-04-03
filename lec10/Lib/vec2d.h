// This is the header file
// of our vec2d class for ROOT

#ifndef ROOT_vec2d
#define ROOT_vec2d

#include "TObject.h"
#include "TBuffer.h"
#include <cmath>

//definition of our vec2d class
class vec2d : public TObject
{
  protected:

  double x;
  double y;

  public:

  // constructor
  vec2d(); //default
  vec2d(double mx, double my);
  // for how to better fill the event class
  // a copy constructor would be good ;-)

  // destructor
  virtual ~vec2d();

  //operators
  vec2d operator += (vec2d &);
  vec2d operator + (vec2d &);

  // Setter functions
  void SetX(double mx) {x=mx;}
  void SetY(double my) {y=my;}

  // Getter functions
  double GetX() {return x;}
  double GetY() {return y;}

  //our member functions
  void Print();

  // small redefition (const retunr)
  // for function used to sort !
  // (make it also inline since it is
  // a simple operation!)

  double Length() const {return sqrt(x*x+y*y);}

  // Functions needed to use Sorting algorithms
  // in ROOT Container Classes!

  Bool_t IsSortable() const {return kTRUE;}

  Int_t Compare(const TObject *obj) const
  {
    if (Length() == ((vec2d*) obj)->Length())
      return 0;
    else if (Length() > ((vec2d*) obj)->Length())
      return -1;
    else if (Length() < ((vec2d*) obj)->Length())
      return 1;
    else
      return -9999;
  }

  ClassDef(vec2d,1); 
};
 
#endif //ROOT_vec2d
