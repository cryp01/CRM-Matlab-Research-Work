function  f  = fCGYM ( r,q,sigma,t,T,C,G,Y,M,w)
%function used in the computation of the coefficient c
%   f is the inverse of the Fourier transformation
%   r : risk free rate
%   q : dividend
%   sigma : volatility
%   t : present time
%   T : maturity
%   w : log value of the asset

s = -C*gamma(-Y)*((M-1)^Y-M^Y+(G+1)^Y-G^Y);
f = exp(-1i*w*(r-q+s)*(T-t))*exp(C*gamma(-Y)*((M+1i*w)^Y-M^Y+(G-1i*w)^Y-G^Y)*(T-t));


end

