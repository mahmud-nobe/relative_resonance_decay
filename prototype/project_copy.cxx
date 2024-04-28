#include <iostream>
#include <fstream>
#include <random>
#include <vector>
#include <TRandom.h>
#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TROOT.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TVector3.h"


// Define the resonance particle and its properties
const double M_RES = 1.51942; // Resonance mass in GeV/c^2
const double WIDTH_RES = 0.01573; // Resonance width in GeV/c^2
//const double LIFETIME_RES = 2.632e-10; // Resonance lifetime in seconds
//const int PDG_RES = 3122; // PDG code of the resonance (Lambda)

// Define the decay products and their properties
const double M_NEUTRON = 0.9395654; // NEUTRON mass in GeV/c^2
const double M_KAON = 0.497611; // KAON mass in GeV/c^2
//const int PDG_NEUTRON = 2212; // PDG code of the NEUTRON
//const int PDG_KAON = 211; // PDG code of the KAON

// Define the collision energy
const double SQRT_SNN = 200.0; // Center-of-mass energy in GeV

// Define the momentum resolution
//const double MOM_RES = 0.0000001; // Momentum resolution (1%)


class FourMomentum {
   private:
      TLorentzVector vec;
      int pdgCode;

   public:
   // Constructor
      FourMomentum(double px, double py, double pz, double e, int pdgCode_) : vec(px, py, pz, e), pdgCode(pdgCode_) {}
      FourMomentum(TLorentzVector v, int pdgCode_) : vec(v), pdgCode(pdgCode_) {}
	
	// Getters
      TLorentzVector L_vec() const {return vec;}
      double Px() const { return vec.Px(); }
      double Py() const { return vec.Py(); }
      double Pz() const { return vec.Pz(); }
      double E() const { return vec.E(); }
      
      // Setters
    void SetPx(double px) { vec.SetPx(px); }
    void SetPy(double py) { vec.SetPy(py); }
    void SetPz(double pz) { vec.SetPz(pz); }
    void SetE(double e) { vec.SetE(e); }
    void SetPdgCode(int pdgCode_) { pdgCode = pdgCode_; }

    // Other functions
    double Pt() const { return vec.Pt(); } // sqrt(vec.Px()^2 + vec.Py()^2 )
    double P() const { return vec.P(); } // sqrt( vec.Px()^2 + vec.Py()^2 + vec.Pz()^2 )
    // dot product = invariant mass
    double M() const { return vec.M(); } // sqrt( vec.E()^2 - vec.P()^2 )
	FourMomentum Boost(const FourMomentum& resonance);
	
	//operators
	FourMomentum operator + (FourMomentum &);
	FourMomentum operator - (FourMomentum &);

      void Print() const {
         std::cout << "(" << Px() << ", " << Py() << ", " << Pz() << ", " << E() << ")" << std::endl;
      }
};

FourMomentum FourMomentum::Boost(const FourMomentum& resonance) 
{
	// get the boost vector from the resonace fourMomentum
    // boost_v = Px/E, Py/E, Pz/E
    TVector3 boost_v = resonance.L_vec().BoostVector();
         
    // Make a copy of the current Lorentz Vector 
    TLorentzVector new_vec = vec;
         
    // Boost this Lorentz Vector using boost_v
	new_vec.Boost(boost_v); // in place boosted momentum
	
    return FourMomentum(new_vec, pdgCode);
}
      

FourMomentum FourMomentum::operator+(FourMomentum &m)
{
	return FourMomentum(vec.Px()+m.Px(),vec.Py()+m.Py(),vec.Pz()+m.Pz(),vec.E()+m.E(), pdgCode);
}

FourMomentum FourMomentum::operator-(FourMomentum &m)
{
	return FourMomentum(vec.Px()-m.Px(),vec.Py()-m.Py(),vec.Pz()-m.Pz(),vec.E()-m.E(), pdgCode);
}



// Function to perform the resonance decay
std::vector<FourMomentum> ResonanceDecay(const FourMomentum& resonance) {
    std::vector<FourMomentum> daughters;

    // Input: the momentum of the resonance in its rest frame
    // Output: the momentum of daughter particles in lab frame
    
    //FourMomentum resonance_lab = resonance.Boost(beta_res, gamma_res);
    double P_res = resonance.P();
    double E_res = resonance.E();
    double M = resonance.M();
   
    // Generate random decay angles in the rest frame
    // rand()/RAND_MAX --> any random number between 0 and 1
    double theta = acos(2.0*rand()/RAND_MAX - 1.0); // arccos( uniform distribution between -1 and 1 )
    double phi = 2.0*3.14159*rand()/RAND_MAX;  // uniform distribution between 0 and 2pi
    double eta = rand()/RAND_MAX - 0.5 // eta [-0.5, 0.5]

    // Calculate the momentum of the daughter particles in the rest frame
    double P_d = sqrt((M*M - M_NEUTRON*M_NEUTRON - M_KAON*M_KAON)*(M*M - M_NEUTRON*M_NEUTRON - M_KAON*M_KAON)
            - 4.0*M_NEUTRON*M_NEUTRON*M_KAON*M_KAON)/(2.0*M_RES);
    double E_d_NEUTRON = sqrt(P_d*P_d + M_NEUTRON*M_NEUTRON);
    double E_d_KAON = sqrt(P_d*P_d + M_KAON*M_KAON);
    
    double P_T = P_d / cosh(eta);

    // Calculate the 4-momenta of the daughter particles in the rest frame
    FourMomentum NEUTRON(P_T*cos(phi), P_T*sin(phi), P_T*sinh(eta), E_d_NEUTRON, PDG_NEUTRON);
    FourMomentum KAON(-P_T*cos(phi), -P_T*sin(phi), -P_T*sinh(eta), E_d_KAON, PDG_KAON);

    //cout<< KAON.Pz()<<" " <<endl;
    // Boost the daughter particles to the lab frame
    FourMomentum NEUTRON_lab = NEUTRON.Boost(resonance);
    FourMomentum KAON_lab = KAON.Boost(resonance);
   //cout<<NEUTRON_lab.Pz()<<endl;
    
    daughters.push_back(NEUTRON_lab);
    daughters.push_back(KAON_lab);

    return daughters;
}


int project_copy(int n_events = 10000, int n_bg = 10) {

	// Define the number of events and background particles
	const int N_EVENTS = n_events; // Number of events	
	const int N_BG = n_bg; // Number of background particles per event

    // Initialize histograms for the invariant mass and background
    TFile *f=new TFile("lorenz.root","RECREATE");
    TH1F* h_mass = new TH1F("h_mass", "Invariant Mass", 100, 1.08, 1.15);
    TH1F* h_bg = new TH1F("h_bg", "Background", 100, 1.08, 1.15);

    // Initialize the random number generator
    TRandom3 rng;
    gRandom->SetSeed (0);

    // Loop over events
    for (int i_event = 0; i_event < N_EVENTS; i_event++) {
    
        // Generate a resonance particle
        double m_res = (float)gRandom->BreitWigner(M_RES, WIDTH_RES);
        // double m_res = breit_wigner(M_RES, WIDTH_RES);
        FourMomentum resonance(0,0, SQRT_SNN/(2),sqrt((SQRT_SNN*SQRT_SNN/4)+m_res*m_res), PDG_RES);

        // Decay the resonance particle
       // std::vector
    	std::vector<FourMomentum> daughters = ResonanceDecay(resonance);

	
	
    // Add Gaussian smearing to the daughter particle momenta (data is too sensitive to handle, so didn't do it)
    	for (auto& daughter : daughters) {
        	//daughter.SetPx(rng.Gaus(daughter.Px(), daughter.Px()*0.00005));
        	//daughter.SetPy(rng.Gaus(daughter.Py(), daughter.Py()*0.00005));
        	//daughter.SetPz(rng.Gaus(daughter.Pz(), daughter.Pz()*0.00005));
        	//daughter.SetE(abs(rng.Gaus(daughter.E(), daughter.E()*0.00005)));
    	}
   
   

    	// Calculate the invariant mass of the daughter particles
   	double inv_mass = (daughters[0] + daughters[1]).M();
	
    	// Fill the invariant mass histogram
    	//h_mass->Fill(inv_mass);
    	h_mass->Fill(rng.Gaus(inv_mass, inv_mass*0.005));

    	// Generate uncorrelated background events
    	int j=0;
    	while(j<N_BG){
    		double bg_inv_mass = rng.Uniform(h_mass->GetXaxis()->GetXmin(), h_mass->GetXaxis()->GetXmax());
    		h_bg->Fill(bg_inv_mass);
    		j++;
    	//h_mass->Fill(bg_inv_mass);
	}
    }
	// Draw the histograms

	//TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
	h_mass->SetLineColor(kRed);
	h_mass->GetXaxis()->SetTitle("Invariant Mass [GeV/c^{2}]");
	h_mass->GetYaxis()->SetTitle("Counts");
	h_mass->Add(h_bg);
	h_mass->Draw();
	
	h_bg->SetLineColor(kBlue);
	h_bg->Draw("SAME");
	//f->Write();
	//f->Close();
	// Print the number of signal and background events in the signal region
	double signal = h_mass->Integral(h_mass->FindBin(M_RES - 3.0*WIDTH_RES), h_mass->FindBin(M_RES + 3.0*WIDTH_RES));
	double bg = h_bg->Integral(h_bg->FindBin(M_RES - 3.0*WIDTH_RES), h_bg->FindBin(M_RES + 3.0*WIDTH_RES));
	std::cout << "Number of signal events: " << signal << std::endl;
	std::cout << "Number of background events: " << bg << std::endl;
	
	return 0;
}


/*
This code should generate events for the decay of a Lambda particle resonance, apply Gaussian smearing to the daughter particle momenta, 
and calculate the invariant mass of the daughter particles. It also includes a Breit-Wigner distribution for generating the resonance mass, 
uncorrelated background events, and histograms for the invariant mass and background. 
Finally, it calculates the number of signal and background events in the signal region and prints them to the console.
*/

