function  f  = fGBM (r,q,sigma,t,T,w)
%function used in the computation of the coefficient c
%   f is the inverse of the Fourier transformation
%   r : risk free rate
%   q : dividend
%   sigma : volatility
%   t : present time
%   T : maturity
%   w : log value of the asset


f = exp(-1i*w*(r-q-sigma*sigma*.5)*(T-t)-sigma*sigma*w*w*(T-t)*.5);



end

