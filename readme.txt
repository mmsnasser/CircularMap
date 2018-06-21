%  CircularMap: A numerical implementation of the circular map in MATLAB
%
%  Copyright Mohamed Nasser 2017
%  Please cite this MATLAB functions as:
%
%  When citing this software please mention the URL of the master repository 
%  (https://github.com/mmsnasser/CircularMap), and the paper
%  M.M.S. Nasser,Fast Computation of the Circular Map, Computational Methods 
%  and Function Theory, 15 (2015) 187-223.
%
%
%  PLEASE note that this toolbox contains the files:
%  zfmm2dpart.m
%  fmm2d_r2012a.mexw32
%  fmm2d_r2012a.mexw64
%  pthreadGC2-w32.dll
%  pthreadGC2-w64.dll
%  From the Toolbox:
%  L. G REENGARD AND Z. G IMBUTAS , FMMLIB2D: A MATLAB toolbox for
%  fast multipole method in two dimensions, Version 1.2, 2012.
%  http://www.cims.nyu.edu/cmcl/fmm2dlib/fmm2dlib.html
%  PLEASE also cite the FMMLIB2D toolbox.

% The function 
%        [zet,zetp,cntd,rad] = circmap(et,etp,alpha,n)
% Compute the circular map from a bounded multiply connected domain G of 
% connectivity m+1 bounded by \Gamma_0,\Gamma_1,...,\Gamma_m where \Gamma_0
% is the exterior boundary onto a bounded multiply connected circular domain 
% \Omega bounded by the circles C_0,C_1,...,C_m where C_0 is the exterior circle
% and C_0 is the unit circle. The function also Compute the circular map for
% unbounded multiply connected circular domain of connectivity m bounded by
% \Gamma_1,...,\Gamma_m onto an unbounded multiply connected domain \Omega 
% bounded by the circles C_1,...,C_m. 
% In particular, the function "circmap.m" compute
% 1. cntd: a vector contains the centers of the circles C_0,C_1,...,C_m for
% bounded G; and the centers of the circles C_1,...,C_m for unbounded G.
% 2. rad: a vector contains the radius of the circles C_0,C_1,...,C_m for
% bounded G; and the radius of the circles C_1,...,C_m for unbounded G.
% 3. zet: the parameterization of the boundary of \Omega.
% 4. zetp: the derivative of the parameterization of the boundary of \Omega
% where
% 1. et: the parameterization of the boundary of G.
% 2. etp: the derivative of the parameterization of the boundary of G
% 3. alpha: for bounded G, alpha is a point in G that will be mapped 
%    onto 0 in \Omega; for unbounded G, alpha=inf and it will be mapped
%    onto inf.
% 4. n: the number of discretization points of each boundary component
% 5. koebetol: the tolerance of Koebe iterative method
% 6. gmrestol: the tolerance of GMRES iterative method
% 7. koebemaxit: the maximum number of iterations allowed for Koebe iterative method
% 8. gmresmaxit:the maximum number of iterations allowed for GMRES iterative method
% 9. iprec: for the accurecy of the FMM (see zfmm2dpart.m).


