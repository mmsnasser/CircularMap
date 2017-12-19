function y = intfft (x)
%%
% Suppose that f(t) is a real 2pi-periodic function on the interval [0,2pi]
% and the 1*n vector x is the values of the function f(t) at the n  
% equidistant points (n must be even)
% t_j=(j-1)*2*pi/n,  j=1,2,...,n.
% The function f(t) must satisfies
% \int_(0)^(2pi)f(t)dt=0  and  \int_(0)^(2pi)cos(nt)f(t)dt=0
% The function 
%           y = intfft (x)
% uses fft to find the trigonometric interpolating polynomial that 
% interpolates the function f(t) at the n points t_1,t_2,...,t_n. Then the 
% function intfft computes the values a periodic antiderivative function 
% g(t) at the above points t_j, i.e., 
% g(t) = \int f(t)dt
% and g(2pi)=g(0). 
% Note that g(t) is only one of the antiderivatives of the function f(t) 
% and hence g(t) is not unique.
%
%% Author: Mohamed M S Nasser, v 1.0, 10 December 2017.
%
% Example
% n     =  2^7;
% t     = [0:2*pi/n:2*pi-2*pi/n].';
% x     =  sin(t).*exp(-cos(t));
% xie   =  exp(-cos(t));
% xi    =  intfft (x);
% error =  xi-xie %must be a constant
%
%   Copyright 2009 by Mohamed M S Nasser
%   Last updated   : 12/04/2017 with MATLAB 2016a
%   Email          : mms.nasser@qu.edu.qa, mms.nasser@gmail.com
%
%   Problems or suggestions? Email me.
%%
n        =  length(x);   
C        =  ifft(x);
%%
E        =  zeros(n,1);
for k=2:n/2
    E(k) = C(k)/(-1i*(k-1));
end
for k=n/2+2:n
    E(k) = C(k)/(-1i*(k-1-n));
end
%%
y        =  real(fft(E));
%%
end