function [rho,cv,cp]=PropH2(p,T)
    R=8.3145/0.002; 
    rho=p/(R*T);
    cp=(6.88+0.000066*T+0.279*10^-6*T^2)/2*4.1868*1000;
    cv=cp-R;
end