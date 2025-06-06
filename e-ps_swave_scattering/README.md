This folder contains the code and input file for e-Ps elastic s-wave scattering.  

The file 'swave.f' computes the stationary and trial phase shifts for e-Ps elastic scaterring for the S partial wave using the Kohn, inverse Kohn, 
and complex Kohn methods for either the singlet or triplet spin states. The only external dependency is to two functions in LAPACK.  

The input file is formatted as follows:  
nonlinear parameter values (alpha, gamma, zeta)  
integration point numbers  (ns, nc)  
omega value  
00 for singlet, 01 for triplet  
number of energies  
energies  
printing indicator

