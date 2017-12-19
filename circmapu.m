function   [zet,zetp,cntd,rad]=circmapb(et,etp,n,koebetol,koebemaxit,gmrestol,gmresmaxit,iprec)
%
% The function 
%        circmapu
% will be called from the function "cirmap.m"
% to compute the circular map from an unbounded multiply connected domain G 
% of connectivity m bounded by \Gamma_1,...,\Gamma_m onto an unbounded 
% multiply connected circular domain \Omega bounded by the circles 
% C_1,...,C_m. 
% In particular, the function "circmap.m" compute
% 1. cntd: a vector contains the centers of the circles C_1,...,C_m.
% 2. rad: a vector contains the radius of the circles C_1,...,C_m.
% 3. zet: the parameterization of the boundary of \Omega.
% 4. zetp: the derivative of the parameterization of the boundary of \Omega.
% where
% koebetol: the tolerance of Koebe iterative method
% gmrestol: the tolerance of GMRES iterative method
% koebemaxit: the maximum number of iterations allowed for Koebe iterative method
% gmresmaxit:the maximum number of iterations allowed for GMRES iterative method
% iprec: for the accurecy of the FMM (see zfmm2dpart.m).

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




%%
m  =  length(et)/n;
%%
zet   =  et; 
zetp  =  etp;
for k=1:m
    cntd(k) = getInteriorPoint(et(1+(k-1)*n :k*n,1), etp(1+(k-1)*n :k*n,1));
end
alp   =  cntd(1);
Error =  1;
no_of_iteration = 0;
while (Error> koebetol)
    %%
    no_of_iteration = no_of_iteration+1
    zetold = zet;
    for (k=1:m)
        etm   = zet((k-1)*n+1:k*n);
        etmp  = zetp((k-1)*n+1:k*n);
        zm    = cntd(k);
        [S,Sp,c]  =  cmu(etm,etmp,zm,n,5,[],gmrestol,gmresmaxit);
        wn    =  exp(-i*S);
        wnpep =  -i.*Sp.*wn;        
        % at the boundary k
        zet((k-1)*n+1:k*n)  = wn;
        zetp((k-1)*n+1:k*n) = wnpep;
        cntd(k) = 0;
        % for the other boundaries
        zet_other  = [zet(1:(k-1)*n);zet(k*n+1:m*n)]; 
        zet_otherp = [zetp(1:(k-1)*n);zetp(k*n+1:m*n)]; 
        Cau   = (zet_other-zm).*(fcau(etm,etmp,wn./(etm-zm),zet_other.',n,c).');
        Cau2  = (fcauep(etm,etmp,wnpep,zet_other.',n,c).');
        Caup  =  zet_otherp.*Cau2;
        zet(1    :(k-1)*n)  = Cau(1:(k-1)*n);
        zet(k*n+1:m*n)      = Cau((k-1)*n+1:(m-1)*n);
        zetp(1    :(k-1)*n) = Caup(1:(k-1)*n);
        zetp(k*n+1:m*n)     = Caup((k-1)*n+1:(m-1)*n);
        for j=1:m
            if(j~=k)
                cntd(j)  = (cntd(j)-zm).*(fcau(etm,etmp,wn./(etm-zm),cntd(j),n,c).');
            end
        end
    end        
    cnst1  = (-1/(n*i)).*sum((zet./(et-alp)).*etp./(et-alp));
    cnst2  = (-1/(n*i)).*sum((zet./cnst1-et).*etp./(et-alp));
    zet   =  zet./cnst1-cnst2;
    zetp  =  zetp./cnst1;
    cntd     =  cntd./cnst1-cnst2;
    %%
    for (j=1:m)
        rad(j,1)     = sum(abs(zet((j-1)*n+1:j*n)-cntd(j)))/n;
    end
    %%
    Error          = norm(zet-zetold,inf)
    %% 
    if (no_of_iteration>=koebemaxit)
        'Error of Koebe iterative method'
        Error
        'Convegence failed after maximum iterations allowed'
        break;
    end
    %%
end
%%


%%
end