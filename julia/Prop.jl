function PropAir(p,T)
    alpha=3.653
    beta_cv=-1.337*10^-3
    gamma=3.294*10^-6
    delta=-1.913*10^-9
    epsilon=0.2763*10^-12
    R=286.987
    rho=p/(R*T)
    cv=R*((alpha+beta_cv*T+gamma*T^2+delta*T^3+epsilon*T^4)-1)
    cp=cv+R
    k=0.009748221+5.27354*10^-10*p+5.5243*10^-5*T
    mu=6.93093*10^-6+2.35465*10^-13*p+3.85177*10^-8*T
    return(rho,cv,cp,k,mu)
end
air=PropAir(101000,298.15)
println(air)

function PropWL(T)
    rho=-3.9837*10^-3*T^2+2.1274*T+7.1666*10^2
    cv=-3.0788*10^-4*T^3+0.30602*T^2-1.0086*10^2*T+1.5207*10^4
    k=-7.0674*10^-6*T^2+5.6853*10^-3*T-0.45602
    mu=-3.7024*10^-9*T^3+3.7182*10^-6*T^2-1.2513*10^-3*T+0.14156
    return(rho,cv,k,mu)
end
wl=PropWL(298.15)
println(wl)

function PropWV(p,T)
    R=8.3145/0.018 
    rho=p/(R*T)
    cp=(7.75+0.0016*T-1.6*10^-7*T^2)/18*4.1868*1000
    cv=cp-R
    return(rho,cv,cp)
end
WV=PropWV(101000,298.15)
println(WV)

function PropH2(p,T)
    R=8.3145/0.002 
    rho=p/(R*T)
    cp=(6.88+0.000066*T+0.279*10^-6*T^2)/2*4.1868*1000
    cv=cp-R
    return(rho,cv,cp)
end
h2=PropH2(101000,298.15)
println(h2)

function PropN2(p,T)
    R=8.3145/0.028
    rho=p/(R*T)
    cp=(6.7619+0.000871*T-1.643*10^-7*T^2)/28*4.1868*1000
    #cp=(6.3+0.001819*T-0.345*10^-6*T^2)/28*4.1868*1000
    cv=cp-R
    return(rho,cv,cp)
end
N2=PropN2(101000,298.15)
println(N2)

function PropO2(p,T)
    R=8.3145/0.032
    rho=p/(R*T)
    cp=(6.9647+1.0971*10^-3*T-2.03*10^-7*T^2)/32*4.1868*1000
    #cp=(7.025+0.00185*T-abs((T-300)*0.075/300))/32*4.1868*1000
    cv=cp-R
    return(rho,cv,cp)
end
O2=PropO2(101000,298.15)
println(O2)


(rho,cv,cp)=PropO2(101000,298.15)
println(rho)