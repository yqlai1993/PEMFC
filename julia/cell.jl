using ModelingToolkit,OrdinaryDiffEq
# using Unitful,RecursiveArrayTools


@parameters t,current,T_ca_in,T_an_in,p_ca_in,p_an_in,p_ca_out,p_an_out,Rhum_ca_in,Rhum_an_in,flow_ca_in,flow_an_in
@variables m_O2_ca(t),m_N2_ca(t),m_water_ca(t),m_H2_an(t),m_water_an(t),T_ca(t),T_an(t),T_cell(t)
@derivatives D'~t


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

function PropH2(p,T)
    R=8.3145/0.002 
    rho=p/(R*T)
    cp=(6.88+0.000066*T+0.279*10^-6*T^2)/2*4.1868*1000
    cv=cp-R
    return(rho,cv,cp)
end

function PropN2(p,T)
    R=8.3145/0.028
    rho=p/(R*T)
    cp=(6.7619+0.000871*T-1.643*10^-7*T^2)/28*4.1868*1000
    #cp=(6.3+0.001819*T-0.345*10^-6*T^2)/28*4.1868*1000
    cv=cp-R
    return(rho,cv,cp)
end

function PropO2(p,T)
    R=8.3145/0.032
    rho=p/(R*T)
    cp=(6.9647+1.0971*10^-3*T-2.03*10^-7*T^2)/32*4.1868*1000
    #cp=(7.025+0.00185*T-abs((T-300)*0.075/300))/32*4.1868*1000
    cv=cp-R
    return(rho,cv,cp)
end

function PropWV(p,T)
    R=8.3145/0.018 
    rho=p/(R*T)
    cp=(7.75+0.0016*T-1.6*10^-7*T^2)/18*4.1868*1000
    cv=cp-R
    return(rho,cv,cp)
end

function PropWL(T)
    rho=-3.9837*10^-3*T^2+2.1274*T+7.1666*10^2
    cv=-3.0788*10^-4*T^3+0.30602*T^2-1.0086*10^2*T+1.5207*10^4
    k=-7.0674*10^-6*T^2+5.6853*10^-3*T-0.45602
    mu=-3.7024*10^-9*T^3+3.7182*10^-6*T^2-1.2513*10^-3*T+0.14156
    return(rho,cv,k,mu)
end

function psat(T) "饱和水蒸汽压"
    T=T-273.15
    # a=-2.1794+0.02953*T-9.1837*10^-5*T^2+1.4454*10^-7*T^3
    a=-2.18+0.0295*T-9.18*10^-5*T^2+1.44*10^-7*T^3
    p=10^(a+5)
end

function latent(T)
    T=T-273.15
    H=(45070-41.9*T+3.44*10^-3*T^2+2.54*10^-6*T^3-8.98*10^-10*T^4)/0.018
end
    
function stateEq(m,R,T,volu) "状态方程"
    p=m*R*T/volu
end

function mvap(m_water,p_sat,volu,R_water,T) "气态和液态质量"
    m_water_max=p_sat*volu/R_water/T
    if m_water<=m_water_max
       m_water_vap=m_water
       m_water_liq=0
     else
       m_water_vap=m_water_max
       m_water_liq=m_water-m_water_max 
    end
    return(m_water_vap,m_water_liq)
end

function mflow(y_O2,p_water,p_air,flow_tot) "阴极质量流量"
    M_O2=32*10^-3
    M_N2=28*10^-3
    y_N2=1-y_O2
    M_air=y_O2*M_O2+y_N2*M_N2
    M_water=18*10^-3

    omega=M_water*p_water/M_air/p_air
    flow_air=flow_tot/(1+omega)
    flow_water=flow_tot-flow_air
    
    x_O2=y_O2*M_O2/M_air
    x_N2=y_N2*M_N2/M_air
    flow_O2=flow_air*x_O2
    flow_N2=flow_air*x_N2

    return F=[flow_air,flow_water,flow_O2,flow_N2]
end

function mflow2(p_water,p_H2,flow_tot) "阳极质量流量"
    M_H2=2*10^-3
    M_water=18*10^-3
    omega=M_water*p_water/M_H2/p_H2
    flow_H2=flow_tot/(1+omega)
    flow_water=flow_tot-flow_H2
    return F=[flow_H2,flow_water]
end

function memW(a)
    Den_mem=2000
    M_mem=1.1
    if 0<a<=1
        content_w=0.043+17.81*a-39.85*a^2+36*a^3
    elseif a<=3
        content_w=14+1.4*(a-1)
        else
        content_w=16.8
    end
    MolarDen_w=Den_mem*content_w/M_mem
    return(content_w,MolarDen_w)
end


#定值    
M_H2=2*10^-3
M_O2=32*10^-3
M_N2=28*10^-3
M_water=18*10^-3
R=8.3145#u"J/(mol*K)"
R_O2=R/M_O2
R_N2=R/M_N2
R_H2=R/M_H2
R_water=R/M_water
F=96485

thick_mem=1.275*10^-4
area_mem=280*10^-4

T_amb=298.15
p_amb=101000
y_O2_amb=0.21
y_N2_amb=1-y_O2_amb
M_air_amb=y_O2_amb*M_O2+y_N2_amb*M_N2
h_amb=17

volu_ca=0.01
k_ca_in=0.36*10^-5
k_ca_out=0.22*10^-5
h_ca=2
volu_an=0.005
k_an_in=0.36*10^-5
k_an_out=0.22*10^-5
h_an=2
dH_react=1.196*10^8

Cp_stack=35
m_stack=10/120
num_cell=35

#阴极入口压力
p_sat_ca_in=psat(T_ca_in)
p_water_ca_in=Rhum_ca_in*p_sat_ca_in
p_air_ca_in=p_ca_in-p_water_ca_in

air_ca_in=PropAir(p_air_ca_in,T_ca_in)
Cp_air_ca_in=air_ca_in[3]
wv_ca_in=PropWV(p_water_ca_in,T_ca_in)
Cp_wv_ca_in=wv_ca_in[3]
wl_ca_in=PropWL(T_ca_in)
Cv_wl_ca_in=wl_ca_in[2]

#阴极内部压力
p_sat_ca=psat(T_ca)
(m_water_vap_ca,m_water_liq_ca)=mvap(m_water_ca,p_sat_ca,volu_ca,R_water,T_ca)
p_O2_ca=stateEq(m_O2_ca,R_O2,T_ca,volu_ca)
p_N2_ca=stateEq(m_N2_ca,R_N2,T_ca,volu_ca)
p_water_ca=stateEq(m_water_vap_ca,R_water,T_ca,volu_ca) 
p_air_ca=p_O2_ca+p_N2_ca
p_ca=p_air_ca+p_water_ca
Rhum_ca=p_water_ca/p_sat_ca

O2_ca=PropO2(p_O2_ca,T_ca)
Cp_O2_ca=O2_ca[3]
N2_ca=PropN2(p_N2_ca,T_ca)
Cp_N2_ca=N2_ca[3]
wv_ca=PropWV(p_water_ca,T_ca)
Cp_wv_ca=wv_ca[3]
wl_ca=PropWL(T_ca)
Cv_wl_ca=wl_ca[2]

#阴极入口流量
flow_ca_in=k_ca_in*(p_ca_in-p_ca)
y_O2_ca_in=y_O2_amb
(flow_air_ca_in,flow_water_ca_in,flow_O2_ca_in,flow_N2_ca_in)=mflow(y_O2_ca_in,p_water_ca_in,p_air_ca_in,flow_ca_in)

#阴极出口流量
flow_ca_out=k_ca_out*(p_ca-p_ca_out)
y_O2_ca=p_O2_ca/p_air_ca
(flow_air_ca_out,flow_water_ca_out,flow_O2_ca_out,flow_N2_ca_out)=mflow(y_O2_ca,p_water_ca,p_air_ca,flow_ca_out)
flow_O2_react=num_cell*M_O2*current/4/F
flow_water_gen=num_cell*M_water*current/2/F

#阳极入口压力
p_sat_an_in=psat(T_an_in)
p_water_an_in=Rhum_an_in*p_sat_an_in
p_H2_an_in=p_an_in-p_water_an_in

H2_an_in=PropH2(p_H2_an_in,T_an_in)
Cp_H2_an_in=H2_an_in[3]
wv_an_in=PropWV(p_water_an_in,T_an_in)
Cp_wv_an_in=wv_an_in[3]
wl_an_in=PropWL(T_an_in)
Cv_wl_an_in=wl_an_in[2]

#阳极内部压力
p_sat_an=psat(T_an)
(m_water_vap_an,m_water_liq_an)=mvap(m_water_an,p_sat_an,volu_an,R_water,T_an)
p_H2_an=stateEq(m_H2_an,R_H2,T_an,volu_an)
p_water_an=stateEq(m_water_vap_an,R_water,T_an,volu_an) 
p_an=p_H2_an+p_water_an
Rhum_an=p_water_an/p_sat_an

H2_an=PropH2(p_H2_an,T_an)
Cp_H2_an=H2_an[3]
wv_an=PropWV(p_water_an,T_an)
Cp_wv_an=wv_an[3]
wl_an=PropWL(T_an)
Cv_wl_an=wl_an[2]

#阳极入口流量
flow_an_in=k_an_in*(p_an_in-p_an)
(flow_H2_an_in,flow_water_an_in)=mflow2(p_water_an_in,p_H2_an_in,flow_an_in)

#阳极出口流量
flow_an_out=k_an_out*(p_an-p_an_out)
(flow_H2_an_out,flow_water_an_out)=mflow2(p_water_an,p_H2_an,flow_an_out)
flow_H2_react=num_cell*M_H2*current/2/F

#膜水合流量
(sigma_ca,c_ca)=memW(Rhum_ca)
(sigma_an,c_an)=memW(Rhum_an)
sigma_ave=(sigma_an+sigma_ca)/2
n_d=0.0029*(sigma_ca+sigma_an)^2/4+0.05*sigma_ave-3.4e-19
D_w=1.25e-6*exp(2416*(1/303-1/T_cell))
currentDen=current/area_mem
N_water_mem=n_d*currentDen/F-D_w*(c_ca-c_an)/thick_mem
flow_water_mem=num_cell*N_water_mem*M_water*area_mem

#电压
volt_nst=1.229-8.5e-4*(T_cell-298.15)+R*T_cell*(log(p_H2_an)+0.5*log(p_O2_ca))/(2*F)
volt_act=0.9514-3.12e-3*T_cell+1.87e-4*T_cell*log(current)-7.4e-5*T_cell*log(p_O2_ca/(5.08e6*exp(-498/T_cell)))
#resistivity_mem=(thick_mem/area_mem)*181.6*(1+0.03*currentDen+0.062*(T_cell/303)^2*currentDen^2.5)/((sigma_an-0.634-3*currentDen)*exp(4.18*(1-303/T_cell)))
resistivity_mem=0.01605-3.5e-5*T_cell+8e-5*current
volt_ohm=current*resistivity_mem
volt_cell=volt_nst-volt_act-volt_ohm
power_cell=volt_cell*current
eff_cell=(flow_H2_react/flow_H2_an_in)*(volt_cell/1.48)

eqs=(D(m_O2_ca)~flow_O2_ca_in-flow_O2_ca_out-flow_O2_react,
    D(m_N2_ca)~flow_N2_ca_in-flow_N2_ca_out,
    D(m_water_ca)~flow_water_ca_in-flow_water_ca_out+flow_water_gen-flow_water_mem,
    D(m_H2_an)~flow_H2_an_in-flow_H2_an_out-flow_H2_react,
    D(m_water_an)~flow_water_an_in-flow_water_an_out-flow_water_mem,
    D(T_ca)~(h_ca*(T_cell-T_ca)+Cp_air_ca_in*flow_air_ca_in*T_ca_in-Cp_O2_ca*(flow_O2_ca_out+flow_O2_react)*T_ca-Cp_N2_ca*flow_N2_ca_out*T_ca)/(Cp_O2_ca*m_O2_ca+Cp_N2_ca*m_N2_ca),
    D(T_an)~(h_an*(T_cell-T_an)+Cp_H2_an_in*flow_H2_an_in*T_an_in-Cp_H2_an*(flow_H2_an_out+flow_H2_react)*T_an)/(Cp_H2_an*m_H2_an),
    D(T_cell)~(h_ca*(T_ca-T_cell)+h_an*(T_an-T_cell)+h_amb*(T_amb-T_cell)+Cp_H2_an*flow_H2_react*(T_an-T_cell)+Cp_O2_ca*flow_O2_react*(T_ca-T_cell)-Cp_wv_ca*flow_water_gen*(T_ca-T_cell)+flow_H2_react*dH_react-num_cell*volt_cell*current)/(Cp_stack*m_stack))
cell=ODESystem(eqs,name=:cell)
var=[m_O2_ca=>0.0375,m_N2_ca=>0.00328,m_water_ca=>0.00019,m_H2_an=>0.00118,m_water_an=>0.00009,T_ca=>308.15,T_an=>308.15,T_cell=>298.15]
par=[current=>20,T_ca_in=>298.15,T_an_in=>298.15,p_ca_in=>505000,p_an_in=>808000,p_ca_out=>202000,p_an_out=>202000,Rhum_ca_in=>0.5,Rhum_an_in=>0]
prob=ODEProblem(cell,var,(0.0,1.0),par;jac=true,sparse=true)
sol=solve(prob,Rodas5())
println(sol[1,:])