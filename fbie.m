function  [mu , h ]  =  fbie(et,etp,A,gam,n,iprec,restart,gmrestol,maxit)
% The function 
%        [mu,h]  =  FBIE(et,etp,A,gam,n,iprec,restart,gmrestol,maxit)
% return the unique solution mu of the integral equation 
%               (I-N)mu=-Mgam 
% and the function 
%                h=[(I-N)gam-Mmu]/2,
% where et is the parameterization of the boundary, etp=et', 
% A=exp(-i\thet)(et-alp) for bounded G and by A=exp(-i\thet) for unbounded
% G, gam is a given function, n is the number of nodes in  each boundary
% component, iprec is the FMM precision flag, restart is the maximum number 
% of GMRES method inner iterations, gmrestol is the tolerance of the GMRES   
% method, and maxit is the maximum number of GMRES method outer iterations
%% Author: Mohamed M S Nasser, v 1.0, 10 December 2017.

%%
a        = [real(et.') ; imag(et.')];
m        =  length(et)/n-1;
b1       = [etp./A].';
[Ub1]    = zfmm2dpart(iprec,(m+1)*n,a,b1,1);
Eone     = (Ub1.pot).';
%%
b(1,1) = 0;
for k=2:n
    b(k,1) = (-1)^(k+1)*(1/n)*cot(pi*(k-1)/n);
end
%%
[mu,~,~,~,~]      = gmres(@(x)fB(x),-fC(gam),restart,gmrestol,maxit);
if( nargout == 2 )
    h   = (fC(mu)-fB(gam))./2;
end
%%
function  hx  = fB (x)
    bx2   = [x.*etp./A].';
    [Ubx2]= zfmm2dpart(iprec,(m+1)*n,a,bx2,1);
    Ex    = (Ubx2.pot).';
    hx    =  2.*x-(2/n).*imag(A.*Eone).*x+(2/n).*imag(A.*Ex);
end
%%
function  hx = fC (x)
    bx    = [x.*etp./A].';
    [Ubx] = zfmm2dpart(iprec,(m+1)*n,a,bx,1);
    Ex    = (Ubx.pot).';
    for k=1:m+1
        hLx(1+(k-1)*n:k*n,1) = ifft(fft(b).*fft(x(1+(k-1)*n:k*n,1)));
    end
    hx    = -(2/n).*real(A.*Ex)+(2/n).*real(A.*Eone).*x+hLx;    
end
%%
end