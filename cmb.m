function  [S,Sp]  =  cmb(et,etp,alp,n,iprec,restart,gmrestol,maxit)
%
% The function 
%        cmb
% Compute the circular map from the bounded simply connected domain inside
% \Gamma_0 onto the unit disk. This is required one time for each iteration
% of Koebe iterative method
%% Author: Mohamed M S Nasser, v 1.0, 10 December 2017.

%%
t        =   (0:2*pi/n:2*pi-2*pi/n).';
%%
A          =  et-alp;
Ap         =  etp;
%%
gam  = -log(abs(et-alp));
muh  = -imag(clog(et-alp));
%%
[phi , h ]  =  fbie(et,etp,A,gam,n,iprec,restart,gmrestol,maxit);
%%
S   = phi-muh;
%%
Sp  =  derfft(S-t)+1;
%%
end