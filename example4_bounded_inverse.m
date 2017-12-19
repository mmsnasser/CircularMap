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
'Example 4: Very thin ellipses'
%%
n         =   2^12
ratio     =   0.2;
t         =   (0:2*pi/n:2*pi-2*pi/n).';
%% 
Mat_coef  =   [0                    2         2
               0.5                  0.99     -0.01
              -0.5                  0.99     -0.01
               ];                
cent      =   Mat_coef(:,1);
radx      =   Mat_coef(:,2);
rady      =   Mat_coef(:,3);
m         =   length(cent)-1
%%
alphain     =    0;
for k=1:m+1
    et(1+(k-1)*n:k*n,1)    =  cent(k)+0.5.*(+radx(k).*cos(t)+i*rady(k).*sin(t));
    etp(1+(k-1)*n:k*n,1)   =          0.5.*(-radx(k).*sin(t)+i*rady(k).*cos(t));
end
%%

%%
[zet,zetp,cntd,rad]=circmap(et,etp,alphain,n);
%%

%%
[xeh,yeh]=meshgrid(-1:0.001:1,[-0.8,-0.6,-0.4,-0.37,-0.34,-0.2,0,0.2,0.34,0.37,0.4,0.6,0.8]);
weh=xeh+i.*yeh;
weho = weh-cntd(1);
weh(abs(weho)>=rad(1))=NaN+i*NaN;
weho = weh-cntd(2);
weh(abs(weho)<=rad(2))=NaN+i*NaN;
weho = weh-cntd(3);
weh(abs(weho)<=rad(3))=NaN+i*NaN;
wehv = weh(abs(weh)>=0).';
%%
[xev,yev]=meshgrid([-1.8:0.2:1.8],-1:0.001:1);
wev=xev+i.*yev;
wevo = wev-cntd(1);
wev(abs(wevo)>=rad(1))=NaN+i*NaN;
wevo = wev-cntd(2);
wev(abs(wevo)<=rad(2))=NaN+i*NaN;
wevo = wev-cntd(3);
wev(abs(wevo)<=rad(3))=NaN+i*NaN;
wevv = wev(abs(wev)>=0).';
%%
figure;
hold on
box on
plot(real(wehv),imag(wehv),'or','MarkerSize',0.5)
plot(real(wevv),imag(wevv),'ob','MarkerSize',0.5)
for k=1:m+1
    c_cr    =  zet((k-1)*n+1:k*n,1); c_cr(n+1)  =  c_cr(1);
    plot(real(c_cr),imag(c_cr),'k')
end
plot(0,0,'ok')
axis equal
axis([-1.05  1.05  -1.05   1.05])
%%

%%
zehv       =  fcau (zet,zetp,et,wehv);
zevv       =  fcau (zet,zetp,et,wevv);
%%
figure;
hold on
box on
plot(real(zehv),imag(zehv),'or','MarkerSize',0.5)
plot(real(zevv),imag(zevv),'ob','MarkerSize',0.5)
for k=1:m+1
    c_cr    =  et((k-1)*n+1:k*n,1); c_cr(n+1)  =  c_cr(1);
    plot(real(c_cr),imag(c_cr),'k')
end
plot(real(alphain),imag(alphain),'ok')
axis equal
axis([-1.05  1.05  -1.05   1.05])
%%