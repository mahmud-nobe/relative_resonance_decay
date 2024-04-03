void vec2dTree()
{
  gSystem->Load("../Lib/libEvent.so");
  
  TFile *f=new TFile("vec2d.root","recreate");
  TTree *t1=new TTree("tvec","a simple Tree with our vec2d class");
  
  Int_t ev;
  vec2d *v=0;

  t1->Branch("v",&v);
  t1->Branch("ev",&ev);
  
  // fill the tree
  for (Int_t i=0; i<10000; i++) {

    v=new vec2d(gRandom->Uniform(-20,20),gRandom->Uniform(-20,20));
    ev = i;
    t1->Fill();
    
    delete v;
   }
   // save the Tree heade; the file will be automatically closed
   // when going out of the function scope
   t1->Write();

   f->Write();
   f->Close();
}
