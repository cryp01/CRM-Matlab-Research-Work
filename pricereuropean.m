%%pricer of the european call option
tic;

%first case GBM model

%parameter : are those defined in the page 16
k1=-1;
k2=2;
m=0;
J=4;
JB=3;
S0=100;
K=120;
x=log(S0/K);
r=0.1;
q=0.0;
sigma=0.25;
T=50;
t=0.0;

tmp=0.0;
for k = k1:k2
    %computation of Si*
    Vmk=0.0;
    if(k2<=0)
        Vmk=0.0;
    else
        for j=1:2^(JB-1)
            k1b=max(k1,0.0);
            k2b=max(k2,0.0);
            Vmk = Vmk + I1(k1b/2^m,k2/2^m,k,m,j,JB)-I2(k1b/2^m,k2/2^m,k,m,j,JB);
        end
        Vmk=Vmk*K*2^(m/2)/2^(JB-1);
    end
        
    
    %compuation of the coefficient c(m,k)
    cmk=0.0;
    for j=0:2^(J-1)
        cmk = cmk+fGBM(r,q,sigma,t,T,(2*j+1)*pi*2^m/2^J)*exp(2*pi*1i*k*j/2^J);       %formula given by the FFT
    end
    cmk = 2^(m/2)*real(exp(1i*k*pi/2^J)*cmk)/2^(J-1);
    
    tmp = tmp+cmk*Vmk;
    
end

v = exp(-r*(T-t))*tmp;
temps = toc;
fprintf('for a digital option (GBM case) with m=%d, K=%d, k1=%d, k2=%d, j=%d, jB=%d we have : %.16f\n',m,K,k1,k2,J,JB,v);
BS = S0*normcdf((log(S0/K)+(r-q+sigma*sigma*.5)*(T-t))/(sigma*sqrt(T-t)),0,1)-K*exp(-r*(T-t))*normcdf((log(S0/K)+(r-q-sigma*sigma*.5)*(T-t))/(sigma*sqrt(T-t)),0,1); %black scholes value of the european option
fprintf('Black scholes  european call option : %.16f\n',BS); 
fprintf('error is : %e\n',abs(v-BS));
fprintf('CPU time (seconds) : %f\n',temps);
clear;

