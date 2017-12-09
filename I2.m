function Ival = I2(a,b ,k,m,j,J)
%function I2 used in the formula of the european option p 12
%   input parameter are : a, b two integer as deined p 12
%   k is the coefficient of cnk
%   m is the paremeter of the walvet computation
%   j is the index of the sum in the formula of Vnk
%   J is the integer used in the coefficient


Cj = (2*j-1)*pi/2^J;

Ival = (sin(Cj*(2^m*b-k))-sin(Cj*(2^m*a-k)))/(Cj*2^m);

end

