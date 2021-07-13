function [rho,cv,cp]=PropWV(p,T)
    R=8.3145/0.018; 
    rho=p/(R*T);
    cp=(7.75+0.0016*T-1.6*10^-7*T^2)/18*4.1868*1000;
    cv=cp-R;
    
end