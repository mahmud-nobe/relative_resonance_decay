// This is the header file
// of our event class for ROOT

#ifndef ROOT_event
#define ROOT_event

#include "vec2d.h"

#include "TObject.h"
#include "TBuffer.h"
#include "TObjArray.h"

//definition of our event class
class event : public TObject
{
  protected:

  int eventNumber;
  // our vec2d object array
  TObjArray *vec2dArray;

  public:

  // constructor
  event(); //default
  event(int nEv);

  // destructor
  virtual ~event();

  // Setter Functions
  void SetEventNumber(int nEv) {eventNumber=nEv;}

  // now we want to add a vec2d object ...
  // there are different ways and different issues
  // concerning object ownership (chap 8 and 16)
 
  void AddVec2d(vec2d *mv) {vec2dArray->AddLast(mv);}
  void AddVec2d(double mx, double my);

  // Getter functions
  int GetEventNumber() {return eventNumber;}
  int GetNumberOfTracks() {return (int) vec2dArray->GetEntries();}
  // nicer getter function to avoid accessing the array itself
  vec2d* GetTrack(int nT) {return (vec2d*) vec2dArray->At(nT);}

  //our member functions
  void Print();
  void PrintTracks();

  void DoSort() {vec2dArray->Sort();}

  ClassDef(event,1); 
};
 
#endif //ROOT_vec2d
