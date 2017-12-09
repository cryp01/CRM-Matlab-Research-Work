function Ival = I1( a,b ,k,m,j,J)
%function I1 used in the formula of the european option p 12
%   input parameter are : a, b two integer as deined p 12
%   k is the coefficient of cnk
%   m is the paremeter of the walvet computation
%   j is the index of the sum in the formula of Vnk
%   J is the integer used in the coefficient


Cj = ((2*j-1)/2^J)*pi;
Ival1 = exp(b)*sin(Cj*((b*2^m)-k))-exp(a)*sin(Cj*((a*2^m)-k));
Ival2 = exp(b)*cos(Cj*((b*2^m)-k))-exp(a)*cos(Cj*((a*2^m)-k));
Ival = Cj*2^m*(Ival1+Ival2/(Cj*2^m))/(1+(Cj*2^m)^2);


end

