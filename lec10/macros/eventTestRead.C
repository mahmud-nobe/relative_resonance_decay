void eventTestRead()
{
  gSystem->Load("../Lib/libEvent.so");
  
  TFile *f=new TFile("eventTest.root","READ");
  TTree *t1=(TTree*) f->Get("event");
  
  event *ev=0;

  t1->SetBranchAddress("event",&ev); 

  int nEv=t1->GetEntries();

  cout<<" # events in TTree = "<<nEv<<endl;
  cout<<endl;

  TH1D *hNtracks=new TH1D("hNtracks","",60,0,60);hNtracks->Sumw2();
  TH1D *hVecLength=new TH1D("hVecLength","",80,0,40);hVecLength->Sumw2();

  for (Int_t i=0; i<nEv; i++) 
    {
      //get event
      t1->GetEntry(i);

      // print info for first 10 events
      if (i<10)
	ev->PrintTracks();
      
      hNtracks->Fill(ev->GetNumberOfTracks());
      
      // loop over all tracks in the event ...
      for (int k=0;k<(int) ev->GetNumberOfTracks();k++)
	{
	  hVecLength->Fill(ev->GetTrack(k)->Length());
	}
    }
  
  // normailze to prob. density
  hNtracks->Scale(1/(double) nEv * 1/(double) hNtracks->GetBinWidth(1));
  
  TCanvas *c1=new TCanvas("c1","#1",800,600);

  // fit with possion to check input !!!
  TF1 *fp=new TF1("fp","TMath::PoissonI(x,[0])",0,60);
  fp->SetParameter(0,20);
  
  hNtracks->Fit(fp,"R");

  hNtracks->DrawCopy();
  
  TCanvas *c2=new TCanvas("c2","#2",800,600);
  hVecLength->DrawCopy();

  f->Close();
}
