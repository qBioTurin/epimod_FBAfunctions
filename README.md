
**epimod_FBAfunctions**  is an R package for supporting the general modeling framework [GreatMod](https://qbioturin.github.io/epimod/) for the integration of metabolic models, represented as a Flux Balance Analysis (FBA) problem, into a dynamic model, represented as a set of ordinary differential equations (ODEs).

# GreatMod and FBA

# Functions

1. **FBA4Greatmod.generation** generates the *FBA_greatmod* class in which all the elements needed for defining the FBA model to pass to GreatMod are saved. It takes in input a Matlab file (for instance downloaded from ....).
2. **setDiet** changes the lower constraint of the FBA model considering a specific diet. 
3. **setConstraints**
