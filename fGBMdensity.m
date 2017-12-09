function  f  = fGBMdensity(r,q,sigma,t,T,m,k1,k2,J,x)
%returns the approximation of level m of the density function of the GBM
%dynamic
%the computation of the sum is done from k1 to k2 (k1<k2)
%r,q,sigma,t,T are the parameters of the GBM dynamics
%m,k1,k2 and J are the parmaetersw  used in the computation of Pm(f) p 5
%formula (8)

f=0.0;
for k=k1:k2
     cmk=0.0;
    
    for jp=0:(2^(J-1))
        cmk = cmk+fGBM(r,q,sigma,t,T,(2*(jp+1))*pi*2^m/2^J)*exp(2*pi*1i*k*jp/2^J);
    end
    cmk = 2^(m/2)*real(exp(1i*k*pi/2^J)*cmk)/2^(J-1);
    
    phimk=0.0;
    
    for jp = 1:2^J-1
       phimk = phimk + cos((2*jp+1)*pi*(2^m*x-k)/2^J);     
    end
    
    phimk = 2^(m/2.)*phimk/2^(J-1);
    
    f=f+cmk*phimk;
end



end

