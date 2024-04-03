#include <iostream>
#include "TROOT.h" //always needed
#include "TStyle.h" //usually needed
#include "TH2.h"
#include "TTree.h"
#include "TFile.h"
#include "TRandom.h"

#include "../../Lib/vec2d.h"
#include "../../Lib/event.h"


void eventTest()
{
  TFile *f=new TFile("eventTest.root","recreate");
  TTree *t1=new TTree("event","a simple Tree with our event class");
  
  event *ev=0;

  t1->Branch("event",&ev); 

  Int_t meanTracksPerEvent=25;

  // nice way to debug root classes in
  // memory. Useful to find memory leaks
  // activate by add or create a .rootrc
  // in home directory and add/change:
  // Root.ObjectStat:         1
  // (turn off if no debug is needed,
  // slows down execution!)
  //gObjectTable->Print();

  // fill the tree
  for (Int_t i=0; i<1000; i++) 
    {
      ev=new event(i);     
      
      // distribute the number of tracks per event
      // according to a poisson!!!
      Int_t nTracksPerEvent=gRandom->Poisson(meanTracksPerEvent);

      for (Int_t k=0;k<nTracksPerEvent;k++)
	{
	  // not the most flexible and elegant, but safest method
	  //ev->AddVec2d(gRandom->Uniform(-20,20),gRandom->Uniform(-20,20));

	  // more generic/cleaner notation wise more obvious
	  // but only works since the TObjArray is the owner
	  // so hence we are not allowed to delete v here
	  // (so here explicitly not every new needs a delete
	  // therefore a bit more dangerous ...)
	  // Remark: Try to implement both to be safe!!!!

	  vec2d *v=new vec2d(gRandom->Uniform(-20,20),gRandom->Uniform(-20,20));
	  ev->AddVec2d(v);

	  // adding the delete will cause a crash, since
	  // delete e; with eevent as the owner, trys
	  // to delete too!!!
	  //delete v;
	}

      // sort them by length!!!
      ev->DoSort();

      t1->Fill();

      //DEBUG:
      //ev->Print();
      if (i<10)
	ev->PrintTracks();

      delete ev;   
    }
  
  //DEBUG:
  //gObjectTable->Print();
  
  t1->Write();
  
  f->Write();
  f->Close();
}

int main()
{
    eventTest();
}

