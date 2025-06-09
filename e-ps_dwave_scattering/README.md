This folder contains the code and input file for D-wave e-Ps elastic scattering.  

Files dwave_master.f and mixed_master.f are identical except for the short-range expansion terms, which are computed in the subroutine SQUARE. mixed_master.f contains the fill short-range expansion with the 'mixed symmetry' terms, while dwave_master.f does not. Both of these files compute singlet and triplet D-wave phase shifts for e-Ps elastic scattering using the Kohn, complex Kohn and inverse Kohn methods.  

The input file, Pd.data, is organized as follows:  
- nonlinear parameters (alpha, gamma, zeta)  
- integration point numbers (ns, nc)  
- omgega value
- 00 for singlet, 01 for triplet
- shielding function pwer (nb)
- number of energies
- energies
- printing parameter
