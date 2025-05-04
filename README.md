IRPCODE2
Computes the first‐Born approximation of the scattering amplitude for a single, perfectly spherical scatterer in a homogeneous medium (following Wu & Aki). Then, using Rose’s backscatter model, it determines the number of particles per unit volume (and corresponding material weight) required to achieve a specified backscatter coefficient range. IRPCODE2 assumes the scatterer material is perfectly compressed and spherical when converting particle count to weight.

Modular components
born_approx_final.py
Inputs: mechanical and acoustic properties of the scatterer and medium, plus test frequency
Output: scattering amplitude (first‐Born approximation)

particule_per_unit_vol_final.py
Inputs: desired scattering amplitude and scatterer material density
Outputs: particles per unit volume and material weight needed for that amplitude
