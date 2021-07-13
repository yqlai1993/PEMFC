function [ rho,cv,cp,k,mu ] = PropAir( p,T )
%PROPAIR 此处显示有关此函数的摘要
%   此处显示详细说明

    alpha=3.653;
    beta_cv=-1.337*10^-3;
    gamma=3.294*10^-6;
    delta=-1.913*10^-9;
    epsilon=0.2763*10^-12;
    R=286.987;
    rho=p/(R*T);
    cv=R*((alpha+beta_cv*T+gamma*T^2+delta*T^3+epsilon*T^4)-1);
    cp=cv+R;
    k=0.009748221+5.27354*10^-10*p+5.5243*10^-5*T;
    mu=6.93093*10^-6+2.35465*10^-13*p+3.85177*10^-8*T;
    
end

