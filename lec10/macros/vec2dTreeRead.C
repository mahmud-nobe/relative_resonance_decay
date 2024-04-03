void vec2dTreeRead()
{
  gSystem->Load("../Lib/libEvent.so");
  
  TFile *f=new TFile("vec2d.root","read");
  TTree *t1=(TTree*) f->Get("tvec");
  
  Int_t ev;
  vec2d *v=0;

  t1->SetBranchAddress("v",&v);
  t1->SetBranchAddress("ev",&ev);
  
  // read the tree
  // only the fist 10 entries
  Int_t nentries = (Int_t)t1->GetEntries();

  cout<<" # entries = "<<nentries<<endl;

  for (Int_t i=0; i<10; i++) {

    t1->GetEntry(i);
    cout<<" Event # = "<<ev<<endl;
    v->Print();

   }
}
