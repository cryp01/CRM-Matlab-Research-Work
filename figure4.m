%script of the plot : figure 4
% we plot the absolute difference between the SWIFT method and the exact
% solution of the cas or nothing options call.
% the aprameters used are : S0 =100, r=0.1, T=0.1, sigma=0.25
% for different strike K (80, 100, 120)


%parameter
S0=100;
r=0.1;
q=0.0;
sigma=0.25;
T=0.1;
t=0.0;





%% we compute the Black scholes value of the cash or nothing call :
K=80;
BS80 = exp(-r*(T-t))*normcdf((log(S0/K)+(r-q-sigma*sigma*.5)*(T-t))/(sigma*sqrt(T-t)),0,1); %black scholes value of the option

K=100;
BS100 = exp(-r*(T-t))*normcdf((log(S0/K)+(r-q-sigma*sigma*.5)*(T-t))/(sigma*sqrt(T-t)),0,1); %black scholes value of the option

K=120;
BS120 = exp(-r*(T-t))*normcdf((log(S0/K)+(r-q-sigma*sigma*.5)*(T-t))/(sigma*sqrt(T-t)),0,1); %black scholes value of the option

x=log(S0/K);
%%we compute the cash or nothing price with the SWIFT method
diff80 = zeros(1,5);    %to store the absolute difference with strike K = 80
diff100 = zeros(1,5);    %to store the absolute difference with strike K = 100
diff120 = zeros(1,5);    %to store the absolute difference with strike K = 120
m=1:5;
k1=[-1;-3;-7;-8;-17];
k2=[2;2;4;16;32];
J = 4:8;
JB=3:7;
for i=1:5
    tmp=0.0;
    for k = k1(i):k2(i)
        %computation of Si*
        Si = 0.0;
        for jp=0:(2^(JB(i)-1)-1)
            Si = Si + sin((jp+.5)*pi*abs(k)/(2^(JB(i)-1)))/(2*jp+1);
        end
        Si=Si*2./pi;

        %computation of V
        Vmk=(sign(k)*Si+.5)/2^(m(i)/2);
        %compuation of the coefficient c(m,*ùk) we use the FFT approach
        cmk=0.0;

        for jp=0:(2^J(i)-1)
            cmk = cmk+fGBM(r,q,sigma,t,T,(2*jp+1)*pi*2^m(i)/2^J(i))*exp(2*pi*1i*k*jp/2^J(i));
        end

        cmk = 2^(m(i)/2)*real(exp(1i*k*pi/2^J(i))*cmk)/2^(J(i)-1);    
        tmp = tmp+cmk*Vmk;    
    end

    v = exp(-r*(T-t))*tmp;
    diff80(i)=abs(v-BS80);
    diff100(i)=abs(v-BS100);
    diff120(i)=abs(v-BS120);
end

figure('Name','absolute errors','NumberTitle','off');
plot(m,diff80,'k',m,diff100,'--k',m,diff120,':k');
xlabel('m');
ylabel('absolute error');
ylim([0 1]);
legend('K=80','K=100','K=120');

%%we plot the absolute difference

