using ModelingToolkit,OrdinaryDiffEq

@parameters t,height1,width1,
@variables T1(t),T2(t),T3(t),m1_air(t),m1_water_vap(t),m2_water(t),m3_hydr(t),m3_water_vap(t)
@derivatives D'~t

function psat(T)
    T=T-273.15
    # a=-2.1794+0.02953*T-9.1837*10^-5*T^2+1.4454*10^-7*T^3
    a=-2.18+0.0295*T-9.18*10^-5*T^2+1.44*10^-7*T^3
    p=10^(a+5)
end

function PropAir(p_air,T_air)
    alpha=3.653
    beta_cv=-1.337*10^-3
    gamma=3.294*10^-6
    delta=-1.913*10^-9
    epsilon=0.2763*10^-12
    R_air=286.987
    rho_air=p_air/(R_air*T_air)
    cv_air=R_air*((alpha+beta_cv*T_air+gamma*T_air^2+delta*T_air^3+epsilon*T_air^4)-1)
    cp_air=cv_air+R_air
    k_air=0.009748221+5.27354*10^-10*p_air+5.5243*10^-5*T_air
    mu_air=6.93093*10^-6+2.35465*10^-13*p_air+3.85177*10^-8*T_air
    return(rho_air,cv_air,cp_air,k_air,mu_air)
end

function Propwater(T_water)
    rho_water=-3.9837*10^-3*T_water^2+2.1274*T_water+7.1666*10^2
    cv_water=-3.0788*10^-4*T_water^3+0.30602*T_water^2-1.0086*10^2*T_water+1.5207*10^4
    k_water=-7.0674*10^-6*T_water^2+5.6853*10^-3*T_water-0.45602
    mu_water=-3.7024*10^-9*T_water^3+3.7182*10^-6*T_water^2-1.2513*10^-3*T_water+0.14156
    return(rho_water,cv_water,k_water,mu_water)
end

function coeffheat(switch,T,p,height,width)
    diam=height*width/(2*(height+width))
    if switch==1
        (Den,Cv,Cp,thcon,dyvis)=PropAir(p,T)
        else
        (Den,Cv,Cp,thcon,dyvis)=PropH2(p,T)
    end
    Re=Den*vel*diam/dyvis
    Pr=Cp*dyvis/thcon
    Nu=0.023*Re^0.8*Pr^0.4 
    h=thcon*Nu/diam 
end

function humidity(switch,T,p,height,width,length,length2)
    if switch==1
        M_O2=32*10^-3
        M_N2=28*10^-3
        x_O2=0.21
        x_N2=1-x_O2
        M_medium=x_O2*M_O2+x_N2*M_N2
        else
        M_medium=2*10^-3
    end
    M_water=18*10^-3
    volu_medium=height*width*length
    volu_water_vap=height*width*length2
    coeff_diff=1.55*T^1.5*sqrt(1/M_medium+1/M_water)/(p*(volu_medium^(1/3)+volu_water_vap^(1/3))^2)
    R_mem=
    thick_mem=
    p_water_vap_cen=
    p_water_vap_mem=
    n_water_vap=coeff_diff*p*(p_water_vap_cen-p_water_vap_mem)/(R_mem*T*p_medium*thick_mem)
    area=height*width
    flow_water_vap_in=n_water_vap*area
    hum_out=flow_water_vap_out/flow_mudiem_out
end

function latentheat(T)
    laheat=-1.2e-3*T^2-1.6485*T+3040.3
end

#控制体1湿空气通道
h1=coeffheat("air",T1,p1,height,width)
area=height*width
Q1=h1*area*(T2-T1)
hum1=humidity(1)
p1_water_vap=hum1*p1/(0.622+hum1)
p1_sat=psat(T1)
Rhum1=p1_water_vap/p1_sat

#控制体3湿氢气通道
h3=coeffheat("H2",T3,p3,height3,width3)
Q3=h3*area*(T2-T3)
hum3=humidity(3)
p3_water_vap=hum3*p3/(8.938+hum3)
p3_sat=psat(T3)
Rhum3=p3_water_vap/p3_sat

#控制体2水通道
latent1=latentheat(T1)
latent3=latentheat(T3)


eqs=(D(T1)~(Q1+mflow1_air_in*H1_air_in+mflow1_water_vap_in*H1_water_vap_in-mflow1_air_out*H1_air_out-mflow1_water_vap_out*H1_water_vap_out)/(C_air*m1_air+C_water_vap*m1_water_vap),
    D(m1_air)~mflow1_air_in-mflow1_air_out,
    D(m1_water_vap)~mflow1_water_vap_in-mflow1_water_vap_out,
    D(T3)~(Q3+mflow3_hydr_in*H3_hydr_in+mflow3_water_vap_in*H3_water_vap_in-mflow3_hydr_out*H3_hydr_out-mflow3_water_vap_out*H3_water_vap_out)/(C_hydr*m3_hydr+C_water_vap*m3_water_vap),
    D(m3_hydr)~mflow3_hydr_in-mflow3_hydr_out,
    D(m3_water_vap)~mflow3_water_vap_in-mflow3_water_vap_out,
    D(T2)~(mflow2_water_in*H2_water_in-mflow2_water_out*H2_water_out-mflow1_water_vap_in*latent1-mflow3_water_vap_in*latent3-Q1-Q3)/(C_mem*m_mem+Cp_water*m2_water),
    D(m2_water)~mflow2_water_in-mflow2_water_out-mflow1_water_vap_in-mflow3_water_vap_in)
varhum=ODESystem(eqs,name=:varhum)