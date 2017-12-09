%script to obtain figure 3
%we plot the density function and modulus of the characteristic function
%corresponding to the CGMY and GBM dynamics
%for GBM the parameter is m=2
% parameters : S0=100, K=100, r=0.1, T=1, sigma=0.25

%plot of the characterisctic function corresponding to the CGMY
%dynamics

%parameters
S0=100;
K=100;
r=0.1;
T=1;
sigma=0.25;
q=0.0;
t=0.0;
C=1.;
G=5.;
M=5;
Y=1.5;
Y1=0.1;



%x axis
x=-10:10;
y=zeros(1,size(x,2));
z=zeros(1,size(x,2));
w=zeros(1,size(x,2));
for i=1:size(x,2)
    y(i)=abs(fGBM(r,q,sigma,t,T,x(i)));
    z(i)=abs(fCGYM(r,q,sigma,t,T,C,G,Y,M,x(i)));
    w(i)=abs(fCGYM(r,q,sigma,t,T,C,G,Y1,M,x(i)));
end
figure('Name','modulus of the characteristic function','NumberTitle','off');
plot(x,y ,':r',x,z,'k',x,w,'--g');
xlabel('x');
ylabel('Fourier transform');
legend('GBM m =2','CGYM Y=1.5 m=0','CGYM Y=.1 m=4');
hold off;

%plot of the density function corresponding to the GBM and the CGYM
%dynamics

x=-2:.1:2;
y=zeros(1,size(x,2));
z=zeros(1,size(x,2));
w=zeros(1,size(x,2));
for i=1:size(x,2)
    y(i)=fGBMdensity(r,q,sigma,t,T,2,-3,2,5,x(i));
    z(i)=fCGYMdensity(r,q,sigma,t,T,C,G,Y,M,0,-1,2,4,x(i));
    w(i)=fCGYMdensity(r,q,sigma,t,T,C,G,Y1,M,4,-8,16,7,x(i));
end
figure('Name','density function','NumberTitle','off');
plot(x,y ,':r',x,z,'k',x,w,'--g');
xlabel('x');
ylabel('density');
legend('GBM m =2','CGYM Y=1.5 m=0','CGYM Y=.1 m=4');
hold off;