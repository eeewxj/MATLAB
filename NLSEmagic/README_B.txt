NLSEmagic1D Binaries Readme
----------------------------
NLSE1D_TAKE_STEPS_XX_XX.cX
Program to integrate a chunk of time steps of the 1D Nonlinear Shrodinger Equation
i*Ut + a*Uxx - V(x)*U + s*|U|^2*U = 0

Ronald M Caplan
Computational Science Research Center
San Diego State University

INPUT:
(U,V,s,a,h2,BC,chunk_size,k)
U  = Current solution matrix
V  = External Potential matrix
s  = Nonlinearity paramater
a  = Laplacian paramater
h2 = Spatial step size squared (h^2)
BC = Boundary condition selection switch:  1: Dirchilet 2: MSD 3: Lap=0 4: One-sided diff
chunk_size = Number of time steps to take
k  = Time step size

OUTPUT:
U:  New solution matrix

Typical call:
U = NLSE1D_TAKE_STEPS_CD(U,V,s,a,h^2,BC,chunk_size,k);

-------------------------------------

Basic installation instructions (for Windows and Linux):
1) Extract all files in this zip file into a directory of your choice.
2) Run MATLAB and change the directory to where you unzipped the files.
3) Now you can run NLSEmagic1D.

Alter the m file NLSEmagic1D.m to run your own simulations.

For examples on how to plot, save movies/images, run timings, etc. 
see the full research script code package of NLSE1D at www.nlsemagic.com.

-------------------------------------------------------------------------
Disclaimer - This code is given as "as is" and is not guaranteed at all 
             and I take no liability for any damage it may cause. Feel 
			 free to distribute the code as long as you keep the authors 
			 name in it.  If you use these codes, the author would appreciate 
			 your acknowledging having done so in any reports, publications, etc. 
			 resulting from their use.  
			 Also, donations are welcome at www.nlsemagic.com.
