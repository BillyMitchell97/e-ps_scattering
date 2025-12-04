This folder contains the code and input files to compute e-Ps phase shifts.  

pwave.f is code for singlet and triplet P-wave phase shifts for e-Ps elastic scattering using the Kohn, inverse Kohn, complex Kohn S-matrix, or complex Kohn T-matrix methods. The only external dependencies are to LAPACK routines.   

The input file, Pw.data, is structured as follows:  
nonlinear parameters (alpha, gamma, zeta)  
integration points (ns, nc)  
omega value  
00 for singlet, 01 for triplet  
shielding function power (nb)  
number of energies  
energies  
print indicator

Example compile command for the Intel IFORT Fortran compiler:
ifort pwave.f -o PWAVE -mkl -O0 -traceback
