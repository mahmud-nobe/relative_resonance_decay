# Final Project
## Relativec Resonance Decay
Final Project for Computational Physics Course to measure the resonance decay of a particle generator.

### Group 2 Members: 
- FNU Syed Faraz Mohd Mehdi
- Robert Dolan
- Wooyong Han
- Mahmud Un Nobe
  
### Relative Resonance Decay:

#### Kinematic information:

pT = 1,2,3 GeV/c of the resonance (or just a uniform distribution, for test purposes, I would start with a fixed pT
mass of resonance and decay particles according to PDG book. Assume both decay into a proton and a π−.
Lab frame: eta=[-0.5,0.5] (uniform), phi=isotropic, sqrt(sNN) = 200 GeV/c

(Lambda(1520) resonance)


#### Implement a relativistic resonance decay
**Monte Carlo Simulation**:
- Use the Lambda particle resonance (look up in PDG)
- Implement decay in resonance rest frame (isotropic decay)
- Boost in the center-of-mass frame of the collision using sqrt(sNN)=200 GeV/c
- Save this event (of course lots of them)

- Write an analysis macro: Invariant Mass calculation: to identify/confirm the resonance
- Additions:
i) Implement Gaussian smearing of uncertainties in momentum (mimic detector resolution)
ii) Add uncorrelated background

**Hints:**
- We need a 4 momentum vector (look at TLorentzVector),
- beneficial to extend this by creating our own track class allowing to store further information like
particle type
- A good starting point for an event class we did in Lecture 10

**Preparations:**
- What is needed and how do you want to store/generate?
  - Track class, Event class, containers/trees ...?
- Prepare a short outline of what functionalities the classes need and
- a quick flow chart of the project itself

