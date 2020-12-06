%%% 23-11-20        first revision
%%% Benoit Duchet, University of Oxford
*************************************************
This package can be used to compute the isostable and Hilbert amplitude fields 
associated with the basin of attraction of a two-dimensional dynamical system 
with complex conjugate eigenvalues. A two-dimensional Wilson-Cowan (WC) model 
is implemented in the relevent functions in the "modules" folder as an example.
*************************************************
*** isostableAmpField_2D_FP.m computes isostable amplitude fields.

*** hilbertAmpField_2D_FP.m computes Hilbert amplitude fields.

*** Reuse the scripts and functions provided for your own applications. 

*** Please send any comment or suggestion to benoit.duchet@ndcn.ox.ac.uk
*************************************************
Note: Non-Windows users should compile ("mex") the c-file provided prior to running the 
example scripts. To this end,
	- set 'modules' as the current folder in Matlab
	- run the following command
			mex fwdSimEulerMaruyama_WC.c

 
			