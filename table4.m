%table 4, we compute the error of the recovered density by the swift method
%at scale m=0 to the CGYM dynamics with parameters S0=100, K=110, r=0.1,
%q=0.05 T=5, C=1, G=5, M=5,Y=1.5


%parameters
S0=100;
K=110;
r=0.1;
q=0.05;
T=5;
t=0.;
C=1;
G=5;
M=5;
Y=1.5;
m=0;

%various interval : 
a=[-1 -2 -5 -10 -20 -32.83];
b=[1 2 5 10 20 25.19];

%we use the formula given table 2 to obtain k1 and k2
k1=[-1 -3 -9 -19 -39 -35]
k2=[2 4 10 20 40 52]

%we compute f(k1/2^m) and f(k2/2^m)

for i=1:size(a,2)
    fprintf('[%f,%f] f(k1/2^m)=%f , f(k2/2^m)=%f and error = %f\n',a(i),b(i),fCGYMdensity(r,q,0,t,T,C,G,Y,M,m,k1(i),k2(i),2*k2(i),k1(i)/2^m),fCGYMdensity(r,q,0,t,T,C,G,Y,M,m,k1(i),k2(i),2*k2(i),k2(i)/2^m),abs(fCGYMdensity(r,q,0,t,T,C,G,Y,M,m,k1(i),k2(i),2*k2(i),k1(i)/2^m)-fCGYMdensity(r,q,0,t,T,C,G,Y,M,m,k1(i),k2(i),2*k2(i),k2(i)/2^m)));
end