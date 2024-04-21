General Steps for Simulations:

1. Choose the mass and width of the resonance particle (res_mass, res_with)

2. Generate the invariant mass of the resonance particle following a breit-Wigner Distribution
	- use gRandom->BreitWigner function from root
	
3. Choose the FourMomentum of the resonance particle in the particle rest frame
	- all momentum are in z axis
	- one choice for $p_z$ = sqrt(sNN)/2
	
4. Assume isotropic decay:
	- eta: from an uniform distribution between -0.5,0.5 (use rand() function)
	- phi: isotropic - uniform distribution from 0 to pi (use rand() function)
	
5. Generate the daughter particles' FourMomentum in the resonance particle rest frame

6. Boost the daughter particles' FourMomentum in the lab frame
	- (physics) get the beta and gamma for the from FourMomentum of resonance particle
	- (code) use the Boost() and BoostVector() function in root to boost the daughter particles' FourMomentum for a given resonace particle's FourMomentum
	
7. Calculate the invariant mass of the daughter particles
	- add their FourMomentum, then take the dot product (Mc^2 = M for c=1)

8. Do step 2-7 for n_events number of time and save the results in a ROOT file and show in the histogram

9. (Optional 1) For each of the event, add an uncorrelated background event
	- randomly select a invariant mass from an uniform distribution around the 
	
