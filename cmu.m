function  [S,Sp,c]  =  cmu(et,etp,alp,n,iprec,restart,gmrestol,maxit)
%
% The function 
%        cmu
% Compute the circular map from the unbounded simply connected domain 
% outside \Gamma_j (for j=1,2,...,m) onto the exterior unit disk. This is 
% required one time (for j=1,2,...,m) for each iteration of Koebe 
% iterative method
%% Author: Mohamed M S Nasser, v 1.0, 10 December 2017.

%%
t        =   (0:2*pi/n:2*pi-2*pi/n).';
%%
A          =  ones(n,1);
Ap         =  zeros(n,1);
%%
gam  =  log(abs(et-alp));
muh  =  imag(clog(et-alp));
%%
[phi , h ]  =  fbie(et,etp,A,gam,n,iprec,restart,gmrestol,maxit);
%%
S   = phi-muh;
ho  = sum(h)/n;
%%
Sp  =  derfft(S-t)+1;
%%
c         =  exp(+ho);
%%
end