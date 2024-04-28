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

///////////////
// Constants //
///////////////

// Define the resonance particle and its properties
const double M_RES = 1.51942; // Resonance mass in GeV/c^2
const double WIDTH_RES = 0.01573; // Resonance width in GeV/c^2
//const double LIFETIME_RES = 2.632e-10; // Resonance lifetime in seconds
const int PDG_RES = 3122; // PDG code of the resonance (Lambda)

// Define the decay products and their properties
const double M_NEUTRON = 0.9395654; // NEUTRON mass in GeV/c^2
const double M_KAON = 0.497611; // KAON mass in GeV/c^2
const int PDG_NEUTRON = 2212; // PDG code of the NEUTRON
const int PDG_KAON = 211; // PDG code of the KAON

// Define the collision energy
const double SQRT_SNN = 200.0; // Center-of-mass energy in GeV

// Define the momentum resolution
//const double MOM_RES = 0.0000001; // Momentum resolution (1%)
const double MOM_WIDTH = 0.005; // std for gaussian smearing of daughter momentum

///////////////////////////////////
// FourMomentum Class Definition //
///////////////////////////////////

class FourMomentum {
   private:
      TLorentzVector vec;
      int pdgCode;

   public:
      // Constructor
      FourMomentum(); // Default Constuctor
      // Constructor with four components and particle type
 	  	FourMomentum(double px, double py, double pz, double e, int pdgCode_) : vec(px, py, pz, e), pdgCode(pdgCode_) {}
 	  	// Constructor with Lorentz vector and particle type
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
    	
    	//operators
			FourMomentum operator + (FourMomentum &);
			FourMomentum operator - (FourMomentum &);

    	// Other functions
    	double Pt() const { return vec.Pt(); } // sqrt(vec.Px()^2 + vec.Py()^2 )
    	double P() const { return vec.P(); } // sqrt( vec.Px()^2 + vec.Py()^2 + vec.Pz()^2 )
    	// dot product = invariant mass
    	double M() const { return vec.M(); } // sqrt( vec.E()^2 - vec.P()^2 )
			FourMomentum Boost(const FourMomentum& resonance);

      void Print() const {
         std::cout << "(" << Px() << ", " << Py() << ", " << Pz() << ", " << E() << ")" << std::endl;
      }
};

///////////////////////////////////////
// FourMomentum Class Implementation //
///////////////////////////////////////

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



/////////////////////////////
// ResonanceDecay Function //
/////////////////////////////

// Input: the FourMomentum of the resonance in its rest frame
// Output: the FourMomentums of daughter particles in lab frame
    
std::vector<FourMomentum> ResonanceDecay(const FourMomentum& resonance, int boost = 1) {
    std::vector<FourMomentum> daughters;
    
    //FourMomentum resonance_lab = resonance.Boost(beta_res, gamma_res);
    double P_res = resonance.P();
    double E_res = resonance.E();
    double M = resonance.M();
   
    // Generate random decay angles in the rest frame
    // rand()/RAND_MAX --> any random number between 0 and 1
    double theta = acos(2.0*rand()/RAND_MAX - 1.0); // cos theta between [-1, 1]
    double phi = 2.0*3.14159*rand()/RAND_MAX;  // uniform distribution between 0 and 2pi

    // Calculate the momentum of the daughter particles in the rest frame
    double P_d = sqrt((M*M - M_NEUTRON*M_NEUTRON - M_KAON*M_KAON)*(M*M - M_NEUTRON*M_NEUTRON - M_KAON*M_KAON) - 4.0*M_NEUTRON*M_NEUTRON*M_KAON*M_KAON)/(2.0*M_RES);
    
    // double P_d = sqrt( ( pow(M_RES, 2) - pow(M_NEUTRON - M_KAON, 2) ) * ( pow(M_RES, 2) - pow(M_NEUTRON + M_KAON, 2) ) ) / 2*M_RES;
    
    double E_d_NEUTRON = sqrt(P_d*P_d + M_NEUTRON*M_NEUTRON);
    double E_d_KAON = sqrt(P_d*P_d + M_KAON*M_KAON);
    

    // Calculate the 4-momenta of the daughter particles in the rest frame
    FourMomentum neutron(P_d*sin(theta)*cos(phi), P_d*sin(theta)*sin(phi), P_d*cos(theta), E_d_NEUTRON, PDG_NEUTRON);
    FourMomentum kaon(-P_d*sin(theta)*cos(phi), -P_d*sin(theta)*sin(phi), -P_d*cos(theta), E_d_KAON, PDG_KAON);

		// Boost the daughter particles to the lab frame
		FourMomentum neutron_lab = neutron.Boost(resonance);
		FourMomentum kaon_lab = kaon.Boost(resonance);

    if (boost) {	// boosted momentum in lab frame			
				daughters.push_back(neutron_lab);
    		daughters.push_back(kaon_lab);
		}
		else{  // momentum in rest frame	 
				daughters.push_back(neutron);
    		daughters.push_back(kaon);
		};
    
    
    return daughters;
}

////////////////////////////////
// Parent Generation Function //
////////////////////////////////

FourMomentum GetParentMomentum(double m_res){

		// Generate random decay angles in the rest frame
    // rand()/RAND_MAX --> any random number between 0 and 1
    double phi = 2.0*3.14159*rand()/RAND_MAX;  // uniform distribution between 0 and 2pi
    double eta = rand()/RAND_MAX - 0.5; // eta [-0.5, 0.5]
    
    double P_T = 2; // 2 GeV/c^2
    double m_T = sqrt( P_T*P_T + m_res*m_res);
    double E_res = sqrt( pow(P_T*cosh(eta), 2) + m_res*m_res ); 
    
    // Calculate the 4-momenta of the resonance particles in the rest frame
    FourMomentum resonance(P_T*cos(phi), P_T*sin(phi), P_T*sinh(eta), E_res, PDG_RES);
   
   	return resonance;
}

///////////////////////////
// Main Project Function //
///////////////////////////

// Input: 
// n_events = total number of events (decay of resonance particle)
// n_bg = total number of uncorrelated background event per each event

// Save the invariant mass of daughter particles generated in each event
// Plot the histogram of invariant mass for both signal and background

int project_copy(int n_events = 10000, int n_bg = 10) {

		// Define the number of events and background particles
		const int N_EVENTS = n_events; // Number of events	
		const int N_BG = n_bg; // Number of background particles per event

    // Initialize histograms for the invariant mass and background
    //TFile *f=new TFile("lorenz.root","RECREATE");
    TH1F* h_mass = new TH1F("h_mass", "Invariant Mass", 100, M_RES - 4*WIDTH_RES, M_RES + 4*WIDTH_RES);
    TH1F* h_mass_rest = new TH1F("h_mass", "Invariant Mass", 100, M_RES - 4*WIDTH_RES, M_RES + 4*WIDTH_RES);
    TH1F* h_bg = new TH1F("h_bg", "Background", 100, M_RES - 4*WIDTH_RES, M_RES + 4*WIDTH_RES);

    // Initialize the random number generator
    TRandom3 rng;
    gRandom->SetSeed (0);

    // Loop over events
    for (int i_event = 0; i_event < N_EVENTS; i_event++) {
    
        // Generate a resonance particle
        double m_res = (float)gRandom->BreitWigner(M_RES, WIDTH_RES);
        FourMomentum resonance = GetParentMomentum(m_res);

        // Decay the resonance particle
    		std::vector<FourMomentum> daughters = ResonanceDecay(resonance);
    		
    		std::vector<FourMomentum> daughters_rest = ResonanceDecay(resonance, 0); // boost = 0  to get the rest frame momentum

	
	
    		// Add Gaussian smearing to the daughter particle momenta
		  	for (auto& daughter : daughters) {
		      	daughter.SetPx(rng.Gaus(daughter.Px(), daughter.Px()*MOM_WIDTH));
		      	daughter.SetPy(rng.Gaus(daughter.Py(), daughter.Py()*MOM_WIDTH));
		      	daughter.SetPz(rng.Gaus(daughter.Pz(), daughter.Pz()*MOM_WIDTH));
		      	daughter.SetE(abs(rng.Gaus(daughter.E(), daughter.E()*MOM_WIDTH)));
		  	}
  

		  	// Calculate the invariant mass of the daughter particles
		 		double inv_mass = (daughters[0] + daughters[1]).M();
		 		double inv_mass_rest = (daughters_rest[0] + daughters_rest[1]).M();
		
		  	// Fill the invariant mass histogram
		  	// h_mass->Fill(inv_mass);
		  	h_mass->Fill(rng.Gaus(inv_mass, inv_mass*MOM_WIDTH));
		  	h_mass_rest->Fill(rng.Gaus(inv_mass, inv_mass*MOM_WIDTH));

		  	// Generate uncorrelated background events
		  	int j=0;
		  	while(j<N_BG){
		  		double bg_inv_mass = rng.Uniform(M_RES - 4*WIDTH_RES, M_RES + 4*WIDTH_RES);
		  		
		  		h_bg->Fill(bg_inv_mass);
		  		j++;
				}
    }
		
		// Draw the histograms

		h_mass->SetLineColor(kRed);
		h_mass->GetXaxis()->SetTitle("Invariant Mass [GeV/c^{2}]");
		h_mass->GetYaxis()->SetTitle("Counts");
		h_mass->Add(h_bg);
		h_mass->Draw();
		
		h_mass_rest->SetLineColor(kOrange);
		h_mass_rest->Add(h_bg);
		h_mass_rest->Draw("SAME");
		
		h_bg->SetLineColor(kBlue);
		h_bg->Draw("SAME");
		
   	auto legend = new TLegend();
   	legend->AddEntry(h_mass,"Resonance Signal","l");
   	legend->AddEntry(h_mass_rest,"Resonance Signal (Rest Frame)","l");
   	legend->AddEntry(h_bg,"Background","l");
   	legend->Draw("SAME");
	
		// TH1F* h_poisson = new TH1F("h_poisson", "Background", 100, M_RES - 4*WIDTH_RES, M_RES + 4*WIDTH_RES);
		
		// for (int k=0; k<1000; k++){
			// double bg_poisson = rng.Poisson(1.5);
			// std:cout << bg_poisson <<endl;
			// h_poisson->Fill(bg_poisson);
		//} 
		
		// h_poisson->Draw();
		//f->Write();
		//f->Close();
		
		// Print the number of signal and background events in the signal region
		double signal = h_mass->Integral(h_mass->FindBin(M_RES - 3.0*WIDTH_RES), h_mass->FindBin(M_RES + 3.0*WIDTH_RES));
		double bg = h_bg->Integral(h_bg->FindBin(M_RES - 3.0*WIDTH_RES), h_bg->FindBin(M_RES + 3.0*WIDTH_RES));
		std::cout << "Number of signal events: " << signal << std::endl;
		std::cout << "Number of background events: " << bg << std::endl;
		std::cout << "Signal to Background Ratio: " << signal/bg << std::endl;
		
		return 0;
}
