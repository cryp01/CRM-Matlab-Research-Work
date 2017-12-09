%script for the figure 6 we plot the graphics of the density function in
%the case of the CGYM dynamics with parameters S0=100, K=110, r=0.1
%q=0.05,T=5,C=1,G=5,M=5,Y=1.5

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

k1=-5;
k2=5;
J=4;

%xaxes

x=-40:40;
y=fCGYMdensity(r,q,0,t,T,C,G,Y,M,m,k1,k2,J,x );

figure('Name','density SWIFT','NumberTitle','off');
plot(x,y);
xlabel('x');
ylabel('density');
legend('SWIFT');
hold off;

xzoom = -10:10;
yzoom=fCGYMdensity(r,q,0,t,T,C,G,Y,M,m,k1,k2,J,xzoom );

figure('Name','density SWIFT zoom','NumberTitle','off');
plot(xzoom,yzoom);
xlabel('x');
ylabel('density');
legend('SWIFT');
hold off;
