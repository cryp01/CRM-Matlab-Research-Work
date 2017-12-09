%second case the CGMY dynmaics

%parameter : are those defined in the page 16
tic;
k1=-1;
k2=2;
m=1;
J=4;
JB=3;
S0=100;
K=80;
x=log(S0/K);
r=0.1;
q=0.0;
sigma=0.25;
T=1;
t=0.0;
C=1;
G=5;
M=5;
Y=0.1;

tmp=0.0;
for k = k1:k2
    %computation of Si*
    Si = 0.0;
    for j=0:2^(JB-1)-1
        Si = Si + sin((j+.5)*pi*abs(k)/(2^(JB-1)))/(2*j+1);
    end
    Si=Si*2./pi;
    %computation of V
    Vmk=(sign(k)*Si+.5)/2^(m/2);
    %compuation of the coefficient c(m,k)
    cmk=0.0;
    for j=1:2^(J-1)
        cmk = cmk+fCGYM(r,q,sigma,t,T,C,G,Y,M,((2j-1)*pi*2^m)/2^J)*exp((1i*k*pi*(2*j-1))/2^j);
    end
    cmk = 2^(m/2)*real(cmk)/2^(J-1);
    
    tmp = tmp+cmk*Vmk;
    
end

v = exp(-r*(T-t))*tmp;
temps = toc;
fprintf('for a digital option (CGYM case) with m=%d, k1=%d, k2=%d, we have : %f\n',m,k1,k2,v);
BS = exp(-r*(T-t))*normcdf((log(S0/K)+(r-q-sigma*sigma*.5)*(T-t))/(sigma*sqrt(T-t)),0,1); %black scholes value of the option
fprintf('Black scholes Digital european option : %f\n',BS); 
fprintf('error is : %f\n',abs(v-BS));
fprintf('CPU time (seconds) : %f\n',temps);
clear;
