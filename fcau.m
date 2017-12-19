function  fz  = fcau(et,etp,f,z,n,finf)
%%
% The function 
%        fz  = FCAU (et,etp,f,z,n,finf)
% return the values of the analytic function f computed using the Cauchy
% integral formula at interior vector of points z, where et is the
% parameterization of the boundary, finf is the values of f at infinity 
% for unbounded G, n is the unber of nodes in each boundary component.
% The integral is discretized using the trapezoidal rule. The summations  
% are computed using the FMM.
%% Author: Mohamed M S Nasser, v 1.0, 10 December 2017.

%%
vz    = [real(z) ; imag(z)];       % target
nz    = length(z);                 % ntarget
a     = [real(et.') ; imag(et.')]; % source
tn    = length(et);                % nsource=(m+1)n
iprec = 5;                         %- FMM precision flag
%%
bf    = [f.*etp].';
[Uf]  = zfmm2dpart(iprec,tn,a,bf,0,0,0,nz,vz,1,0,0);
b1    = [etp].';
[U1]  = zfmm2dpart(iprec,tn,a,b1,0,0,0,nz,vz,1,0,0);
if( nargin == 4 ) 
    fz    = (Uf.pottarg)./(U1.pottarg);
end
%%
if( nargin == 6 ) 
    fz= (finf-(Uf.pottarg)./(n*i))./(1-(U1.pottarg)./(n*i));
end
%%
end