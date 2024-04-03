// Implementation of our event class for ROOT
// now we have to include our header
// file !!!

#include "event.h"

#include <iostream>  

ClassImp(event)

using namespace std;

event::event()
{
  //default constructor
  //DEBUG:
  //cout<<" event: default constructor called"<<endl;

  eventNumber=-1;
  vec2dArray=new TObjArray();

  // Set explicityly that the 
  // array is the owner!
  // This is important to properly
  // allow to delete memory !!!!
  // This is a caveat/inconsitency 
  // in ROOT, so the best way is to
  // play around and try to create
  // an memory leak/error to truly
  // undertsand this aspect!!!
  vec2dArray->SetOwner(kTRUE);

}

event::event(int nEv)
{
  //user defined constructor
  //cout<<" event: user defined constructor called"<<endl;

  eventNumber=nEv;
  vec2dArray=new TObjArray(50);
  vec2dArray->SetOwner(kTRUE);
}

event::~event() 
{
  //cout<<" event: default destructor called"<<endl;
  vec2dArray->Delete();
  delete vec2dArray;
}; 

void event::AddVec2d(double mx, double my)
{
  vec2d *v=new vec2d(mx,my);
  vec2dArray->AddLast(v);
}

void event::Print()
{
  cout<<"Event Class:"<<endl;
  cout<<"Event # = "<<eventNumber<<endl;
  cout<<"# vec2d = "<<vec2dArray->GetEntries()<<endl;
  cout<<endl;
}

void event::PrintTracks()
{
  cout<<"Event Class:"<<endl;
  cout<<"Event # = "<<eventNumber<<endl;
  cout<<"# vec2d = "<<vec2dArray->GetEntries()<<endl; 

  for (int i=0;i<(int) vec2dArray->GetEntries();i++)
    {
      cout<<"Track "<<i<<" with length "<<GetTrack(i)->Length()<<endl;
    }

  cout<<endl;
}
