There are **three** important parts of the codes:

## FourMomentum class: 
To store the four momentum and related functions.

### Members:
1. 4D vector (using TLorentzVector from ROOT)
2. particle type (either particle name as a string, or PDG code as an integer)

### Functions:

#### Comstructors:
- Default Constructor: FourMomentum()
- Constructor with four double numbers for four component and particle code
- Constructir with a 4D Lorentz vector and particle code

#### Setter: 
SetPx, SetPy, SetPz, SetE, SetPDG

#### Getter: 
GetPx, GetPy, GetPz, GetE, GetPDG

#### Operators:
- **+** operator to add two four momentum
- **+=** operator to add two four momentum (in place)
- **-** operator to subtract two four momentum
- **-=** operator to subtract two four momentum (in place)

#### Other Functions:
- Print()
- get\_inv_mass(): dot product of the FourMomentum
- get\_p_tot(): total momentum
- get\_p_transfer(): transfer momentum
- **Boost(FourMomentum &p):** returns the boosted FourMomentum using the lorentz matrix constructed using FourMomentum p.


## Resonace Decay Function: 
- Returns the daughter particles' FourMomentums in lab frame (or maybe invariant mass)
- Input: the FourMoementum of resonance particle in its rest frame
- Steps 4, 5 and 6 from [readme.md](https://github.com/mahmud-nobe/relative_resonance_decay/blob/main/final_project/readme.md)
- save them in an array and return that

## Main Project Function:
- Input: Number of events, Number of background events per signal event

**For Each Event:**
- Generate resonace particle's rest frame FourMomentum
- Get the daughter particles' lab frame FourMomentums using ResonanceDecay function
- calculate the invariant mass
- as a background event, randomly select an invariant mass (or the FourMomentum) from an uniform distribution around our resonance invariant mass (or the FourMoementum)

**Output:**
- Save the results in a ROOT file and show in histogram

