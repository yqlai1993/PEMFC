using ModelingToolkit,OrdinaryDiffEq
using PyCall

@parameters t,p_out,p_in_design,p_out_design,T_in_design,Rspeed_design,flow_design,eff_design,J,eff_motor,kt,kv,volt_motor,Resist_motor
@variables Rspeed(t)
@derivatives D'~t
@pyimport CoolProp.CoolProp as CP

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

T_in_design=298.15
p_in_design=101000
p_out_design=404000
presratio_design=p_out_design/p_in_design
Rspeed_design=5000
flow_design=1.5e-4
eff_design=0.88
J=5*10^-5

T_in=298.15
p_in=101000
Rspeedratio=(Rspeed/sqrt(T_in))/(Rspeed_design/sqrt(T_in_design))
#a^(1/3)>=2/3*b
a=1.06
b=0.36
c1=Rspeedratio/(a*(1-b/Rspeedratio)+Rspeedratio*(Rspeedratio-b)^2)
c2=(a-2*b*Rspeedratio^2)/(a*(1-b/Rspeedratio)+Rspeedratio*(Rspeedratio-b)^2)
c3=-(a*b*Rspeedratio-b^2*Rspeedratio^3)/(a*(1-b/Rspeedratio)+Rspeedratio*(Rspeedratio-b)^2)
c4=0.3
presratio=p_out/p_in
#pi=pi_design*(c1*flowratio^2+c2*flowratio+c3)
flowratio=(-c2+sqrt(c2^2-4*c1*(c3-presratio/presratio_design)))/(2*c1)
#flowratio=(flow*sqrt(T_in)/p_in)/(flow_design*sqrt(T_in_design)/p_in_design)
flow=flowratio*(flow_design*sqrt(T_in_design)/p_in_design)*p_in/sqrt(T_in)
eff=eff_design*((1-c4*(1-Rspeedratio)^2)*(Rspeedratio/flowratio)*(2-Rspeedratio/flowratio))
Cp_in=1004
Cratio=1.4
# (rho_in,Cv_in,Cp_in,k_in,mu_in)=PropAir(p_in,T_in)
# Cratio=Cp_in/Cv_in
angelvel=2*pi*Rspeed/60
MoForce_comp=flow*Cp_in*T_in*(presratio^(1-1/Cratio)-1)/angelvel/eff

eff_motor=0.98
kt=0.0153
kv=0.0153
volt_motor=320
Resist_motor=0.82
MoForce_motor=eff_motor*kt*(volt_motor-kv*angelvel)/Resist_motor

# H_in=CP.PropsSI("H","T",T_in,"P",p_in,"air")
# S_in=CP.PropsSI("S","T",T_in,"P",p_in,"air")
# T_out=T_in+T_in*((presratio)^(1-1/Cratio)-1)/eff
# H_out=CP.PropsSI("H","S",S_in,"T",T_out,"air")
# power_comp=flow*(H_out-H_in)/eff

eqs=(D(Rspeed)~(MoForce_motor-MoForce_comp)/J)
varcomp=ODESystem(eqs,name=:varcomp)

v=[Rspeed=>5000]
p=[p_in=>101000,p_out=>505000,T_in=>298.15]
prob=ODEProblem(varcomp,v,(0,100),p;jac=true,sparse=true)
sol=solve(prob,Rodas5())
println(sol)