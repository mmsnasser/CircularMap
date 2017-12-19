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
clear
clc
%%
n         =   2^10
ratio     =   0.2;
t         =   (0:2*pi/n:2*pi-2*pi/n).';
%% 
Mat_coef  =   [
               0.0                  2.0       0.1
               0.3+0.35i            0.1       0.5
              -0.3+0.35i            0.1       0.5
               0.0+0.20i            0.4       0.1
               0.0+0.40i            0.4       0.1
               0.0+0.65i            0.6       0.1
               ];               
cent      =   Mat_coef(:,1);
radx      =   Mat_coef(:,2);
rady      =   Mat_coef(:,3);
m         =   length(cent)
%%
for k=1:m
    et(1+(k-1)*n:k*n,1)    =  cent(k)+0.5.*(+radx(k).*cos(t)-i*rady(k).*sin(t));
    etp(1+(k-1)*n:k*n,1)   =          0.5.*(-radx(k).*sin(t)-i*rady(k).*cos(t));
end
%%
figure;
hold on
box on
k=1;
for k=1:m
    c_cr    =  et((k-1)*n+1:k*n,1); c_cr(n+1)  =  c_cr(1);
    plot(real(c_cr),imag(c_cr),'b','LineWidth',2.5)
end
c_cr    =  cent(2:end); 
axis equal
axis([-1.05  1.05  -1.05   1.05])
%%

%%
[zet,zetp,cntd,rad]=circmap(et,etp,inf,n);
%%
figure;
hold on
box on
for k=1:m
    c_cr    =  zet((k-1)*n+1:k*n,1); c_cr(n+1)  =  c_cr(1);
    plot(real(c_cr),imag(c_cr),'b','LineWidth',2.5)
end
c_cr    =  cntd(2:end); 
axis equal
%%

