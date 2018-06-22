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
               0.5                  0.8       0.2
               0.5+0.4i             0.3       0.5
               0.2+0.75i            0.6       0.2
               0.1+0.55i            0.4       0.1
              -0.2+0.35i            0.2       0.5
              -0.55+0.5i            0.4       0.3
              -0.80                 0.2       0.7
              -0.4-0.05i            0.5       0.2
              -0.35-0.35i           0.7       0.3
              +0.1-0.35i            0.1       0.5
               0.00-0.8i            0.7       0.2
               0.55-0.4i            0.6       0.2
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
[zet,zetp,cntd,rad]=circmap(et,etp,inf,n);
%%

%%
[xh,yh]=meshgrid(-2:0.001:2,-1.5:0.1:1.5);
zeh=xh+i.*yh;
[mh,nh] =  size(zeh);
for k=1:m
    zho=(zeh-cent(k))./0.5;
    zeh((real(zho)./radx(k)).^2+(imag(zho)./rady(k)).^2<=1)=NaN+i*NaN;
end
z_ind = 1;    
for k=1:mh
    for j=1:nh
        if (abs(zeh(k,j))>=0)
            zehv(z_ind)=zeh(k,j);
            z_ind  = z_ind+1;
        end
    end
end
wehv    = (zehv-cent(1)).*(fcau(et,etp,zet./(et-cent(1)),zehv,n,1));
z_ind    =   1;
for k=1:mh
    for j=1:nh
        if (abs(zeh(k,j))>=0)
            weh(k,j)=wehv(z_ind);
            z_ind  = z_ind+1;
        else
            weh(k,j) = NaN+i*NaN;            
        end
    end
end
%%
[xv,yv]=meshgrid(-2:0.1:2,-1.5:0.001:1.5);
zev=xv+i.*yv;
[mv,nv] =  size(zev);
for k=1:m
    zvo=(zev-cent(k))./0.5;
    zev((real(zvo)./radx(k)).^2+(imag(zvo)./rady(k)).^2<=1)=NaN+i*NaN;
end
z_ind = 1;    
for k=1:mv
    for j=1:nv
        if (abs(zev(k,j))>=0)
            zevv(z_ind)=zev(k,j);
            z_ind  = z_ind+1;
        end
    end
end
wevv    = (zevv-cent(1)).*(fcau(et,etp,zet./(et-cent(1)),zevv,n,1));
z_ind    =   1;
for k=1:mv
    for j=1:nv
        if (abs(zev(k,j))>=0)
            wev(k,j)=wevv(z_ind);
            z_ind  = z_ind+1;
        else
            wev(k,j) = NaN+i*NaN;            
        end
    end
end
%%

%%
figure;
hold on
box on
for k=1:mh
    plot(real(zeh(k,:)),imag(zeh(k,:)),'r')
end
for k=1:nv
    plot(real(zev(:,k)),imag(zev(:,k)),'b')
end
for k=1:m
    c_cr    =  et((k-1)*n+1:k*n,1); c_cr(n+1)  =  c_cr(1);
    plot(real(c_cr),imag(c_cr),'k')
end
axis equal
axis off
set(gca,'LooseInset',get(gca,'TightInset'))
%%
figure;
hold on
box on
for k=1:mh
    plot(real(weh(k,:)),imag(weh(k,:)),'r')
end
for k=1:nv
    plot(real(wev(:,k)),imag(wev(:,k)),'b')
end
for k=1:m
    c_cr    =  zet((k-1)*n+1:k*n,1); c_cr(n+1)  =  c_cr(1);
    plot(real(c_cr),imag(c_cr),'k')
end
axis equal
axis off
set(gca,'LooseInset',get(gca,'TightInset'))
%%


%%
[xhi,yhi]=meshgrid(-2:0.001:2,-1.5:0.1:1.5);
wehi=xhi+i.*yhi;
[mhi,nhi] =  size(wehi);
for k=1:m
    wehi(abs(wehi-cntd(k))<=rad(k)) = NaN+i*NaN;
end
z_ind = 1;    
for k=1:mhi
    for j=1:nhi
        if (abs(wehi(k,j))>=0)
            wehiv(z_ind)=wehi(k,j);
            z_ind  = z_ind+1;
        end
    end
end
zehvi    = (wehiv-cntd(1)).*(fcau(zet,zetp,et./(zet-cntd(1)),wehiv,n,1));
z_ind    =   1;
for k=1:mhi
    for j=1:nhi
        if (abs(wehi(k,j))>=0)
            zehi(k,j)=zehvi(z_ind);
            z_ind  = z_ind+1;
        else
            zehi(k,j) = NaN+i*NaN;            
        end
    end
end
%%
[xvi,yvi]=meshgrid(-2:0.1:2,-1.5:0.001:1.5);
wevi=xvi+i.*yvi;
[mvi,nvi] =  size(wevi);
for k=1:m
    wevi(abs(wevi-cntd(k))<=rad(k)) = NaN+i*NaN;
end
z_ind = 1;    
for k=1:mvi
    for j=1:nvi
        if (abs(wevi(k,j))>=0)
            weviv(z_ind)=wevi(k,j);
            z_ind  = z_ind+1;
        end
    end
end
zevvi    = (weviv-cntd(1)).*(fcau(zet,zetp,et./(zet-cntd(1)),weviv,n,1));
z_ind    =   1;
for k=1:mvi
    for j=1:nvi
        if (abs(wevi(k,j))>=0)
            zevi(k,j)=zevvi(z_ind);
            z_ind  = z_ind+1;
        else
            zevi(k,j) = NaN+i*NaN;            
        end
    end
end
%%

%%
figure;
hold on
box on
for k=1:mhi
    plot(real(wehi(k,:)),imag(wehi(k,:)),'r')
end
for k=1:nvi
    plot(real(wevi(:,k)),imag(wevi(:,k)),'b')
end
for k=1:m
    c_cr    =  zet((k-1)*n+1:k*n,1); c_cr(n+1)  =  c_cr(1);
    plot(real(c_cr),imag(c_cr),'k')
end
axis equal
axis off
set(gca,'LooseInset',get(gca,'TightInset'))
%%
figure;
hold on
box on
for k=1:mhi
    plot(real(zehi(k,:)),imag(zehi(k,:)),'r')
end
for k=1:nvi
    plot(real(zevi(:,k)),imag(zevi(:,k)),'b')
end
for k=1:m
    c_cr    =  et((k-1)*n+1:k*n,1); c_cr(n+1)  =  c_cr(1);
    plot(real(c_cr),imag(c_cr),'k')
end
axis equal
axis off
set(gca,'LooseInset',get(gca,'TightInset'))
%%
