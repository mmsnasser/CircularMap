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
[xh,yh]=meshgrid(-2:0.0001:2,[-0.8,-0.6,-0.4,-0.2,-0.1,-0.05,0,0.05,0.1,0.2,0.4,0.6,0.8]);
zh=xh+i.*yh;
zho = 2.*(zh-cent(1));
zh((real(zho)./radx(1)).^2+(imag(zho)./rady(1)).^2>=1)=NaN+i*NaN;
zho = 2.*(zh-cent(2));
zh((real(zho)./radx(2)).^2+(imag(zho)./rady(2)).^2<=1)=NaN+i*NaN;
zho = 2.*(zh-cent(3));
zh((real(zho)./radx(3)).^2+(imag(zho)./rady(3)).^2<=1)=NaN+i*NaN;
zhv = zh(abs(zh)>=0).';
%%
[xv,yv]=meshgrid([-1.8:0.2:1.8],-1:0.001:1);
zv=xv+i.*yv;
zvo = 2.*(zv-cent(1));
zv((real(zvo)./radx(1)).^2+(imag(zvo)./rady(1)).^2>=1)=NaN+i*NaN;
zvo = 2.*(zv-cent(2));
zv((real(zvo)./radx(2)).^2+(imag(zvo)./rady(2)).^2<=1)=NaN+i*NaN;
zvo = 2.*(zv-cent(3));
zv((real(zvo)./radx(3)).^2+(imag(zvo)./rady(3)).^2<=1)=NaN+i*NaN;
zvv = zv(abs(zv)>=0).';
%%
figure;
hold on
box on
plot(real(zhv),imag(zhv),'or','MarkerSize',0.5)
plot(real(zvv),imag(zvv),'ob','MarkerSize',0.5)
for k=1:m+1
    c_cr    =  et((k-1)*n+1:k*n,1); c_cr(n+1)  =  c_cr(1);
    plot(real(c_cr),imag(c_cr),'k')
end
plot(real(alphain),imag(alphain),'ok')
axis equal
axis([-1.05  1.05  -1.05   1.05])
%%

%%
[zet,zetp,cntd,rad]=circmap(et,etp,alphain,n);
%%
whv       =  fcau (et,etp,zet,zhv);
wvv       =  fcau (et,etp,zet,zvv);
%%
figure;
hold on
box on
plot(real(whv),imag(whv),'or','MarkerSize',0.5)
plot(real(wvv),imag(wvv),'ob','MarkerSize',0.5)
for k=1:m+1
    c_cr    =  zet((k-1)*n+1:k*n,1); c_cr(n+1)  =  c_cr(1);
    plot(real(c_cr),imag(c_cr),'k')
end
plot(0,0,'ok')
axis equal
axis([-1.05  1.05  -1.05   1.05])
%%
