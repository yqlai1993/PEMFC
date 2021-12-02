using ModelingToolkit,OrdinaryDiffEq


@parameters t,p_out,J
@variables Rspeed(t),MoForce_comp(t),MoForce_motor(t)
D=Differential(t)

function Comp(p_out,Rspeed)  
    T_in_design=298.15
    p_in_design=101000
    p_out_design=505000
    presratio_design=p_out_design/p_in_design
    Rspeed_design=5000
    flow_design=1.5e-2
    eff_design=0.88

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
    angelvel=2*pi*Rspeed/60
    MoForce_comp=flow*Cp_in*T_in*(presratio^(1-1/Cratio)-1)/angelvel/eff
    return(MoForce_comp)
end
@register Comp(p_out,Rspeed)

function Motor(Rspeed)
    eff_motor=0.98
    kt=0.0153
    kv=0.0153
    volt_motor=320
    Resist_motor=0.82
    angelvel=2*pi*Rspeed/60
    MoForce_motor=eff_motor*kt*(volt_motor-kv*angelvel)/Resist_motor
    return(MoForce_motor)
end
@register Motor(Rspeed)

# H_in=CP.PropsSI("H","T",T_in,"P",p_in,"air")
# S_in=CP.PropsSI("S","T",T_in,"P",p_in,"air")
# T_out=T_in+T_in*((presratio)^(1-1/Cratio)-1)/eff
# H_out=CP.PropsSI("H","S",S_in,"T",T_out,"air")
# power_comp=flow*(H_out-H_in)/eff

eqs=[MoForce_comp~Comp(p_out,Rspeed),
    MoForce_motor~Motor(Rspeed),
    D(Rspeed)~(MoForce_motor-MoForce_comp)/J]
varcomp=ODESystem(eqs)
varcomp_simple=structural_simplify(varcomp)

v=[Rspeed=>5000]
tspan=(0,100)
p=[p_out=>505000,J=>5*10^-5]
prob=ODEProblem(varcomp_simple,v,tspan,p)
sol=solve(prob,Rodas5())
println(sol)