%%pricer of the digital option
%we fix the niumbers of digital numbers
digits(32);


%first case GBM model

%parameter : are those defined in the page 16
%start time
tic;
m=5;
K=80;
k1=-17;
k2=32;
J=8;
JB=7;
S0=100;
r=0.1;
q=0.0;
sigma=0.25;
T=.1;
t=0.0;
x=log(S0/K);

tmp=0.0;
for k = k1:k2
    %computation of Si*
    Si = 0.0;
    for jp=0:2^(JB-1)-1
        Si = Si + sin((jp+.5)*pi*abs(k)/(2^(JB-1)))/(2*jp+1);
    end
    Si=Si*2./pi;
    %computation of V
    Vmk=(sign(k)*Si+.5)/2^(m/2);
    %compuation of the coefficient c(m,k) we use the FFT approach
    cmk=0.0;
    for jp=0:(2^J-1)
        cmk = cmk+fGBM(r,q,sigma,t,T,(2*jp+1)*pi*2^m/2^J)*exp(2*pi*1i*k*jp/2^J);
        %w= ((2*j+1)*pi*2^m/2^J)
        %exp(-1i*w*(r-q-sigma*sigma*.5)*(T-t)-sigma*sigma*w*w*(T-t)*.5)
    end
    
    cmk = 2^(m/2)*real(exp(1i*k*pi/2^J)*cmk)/2^(J-1);
    
    tmp = tmp+cmk*Vmk;
    
end

v = exp(-r*(T-t))*tmp;
temps = toc;
fprintf('for a digital option (GBM case) with m=%d, K=%d, k1=%d, k2=%d, j=%d, jB=%d we have : %f\n',m,K,k1,k2,J,JB,v);
BS = exp(-r*(T-t))*normcdf((log(S0/K)+(r-q-sigma*sigma*.5)*(T-t))/(sigma*sqrt(T-t)),0,1); %black scholes value of the option
fprintf('Black scholes Digital european option : %f\n',BS); 
fprintf('error is : %f\n',abs(v-BS));
fprintf('CPU time (seconds) : %f\n',temps);
clear;

% %second test
% %parameter : are those defined in the page 16
% %start time
% tic;
% m=2;
% K=120;
% k1=-3;
% k2=2;
% J=5;
% JB=4;
% S0=100;
% r=0.1;
% q=0.0;
% sigma=0.25;
% T=.1;
% t=0.0;
% x=log(S0/K);
% 
% tmp=0.0;
% for k = k1:k2
%     %computation of Si*
%     Si = 0.0;
%     for j=0:2^(JB-1)-1
%         Si = Si + sin((j+.5)*pi*abs(k)/(2^(JB-1)))/(2*j+1);
%     end
%     Si=Si*2./pi;
%     %computation of V
%     Vmk=(sign(k)*Si+.5)/2^(m/2);
%     %compuation of the coefficient c(m,k) we use the FFT approach
%     cmk=0.0;
%     for j=0:2^J-1
%         cmk = cmk+fGBM(r,q,sigma,t,T,(2j+1)*pi*2^m/2^J)*exp(2*pi*1i*k*j/2^J);
%     end
%     cmk = 2^(m/2)*real(exp(1i*k*pi/2^J)*cmk)/2^(J-1);
%     
%     tmp = tmp+cmk*Vmk;
%     cmk
% end
% 
% v = exp(-r*(T-t))*tmp;
% temps = toc;
% fprintf('for a digital option (GBM case) with m=%d, K=%d, k1=%d, k2=%d, j=%d, jB=%d we have : %f\n',m,K,k1,k2,J,JB,v);
% fprintf('error is : %f\n',abs(v-0.5293295436540910));
% fprintf('CPU time (seconds) : %f\n',temps);
% 
