using ModelingToolkit,OrdinaryDiffEq

@parameters t,flow_in,T_out,p_next
@variables m(t),p_out(t)
@derivatives D'~t

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

function pipeout(switch,T_out,p_out,p_next)
    M_O2=32*10^-3
    M_N2=28*10^-3
    x_O2=0.21
    x_N2=1-x_O2
    M_air=x_O2*M_O2+x_N2*M_N2
    R_air=8.3145/M_air
    (rho_air,Cv_air,Cp_air,k_air,mu_air)=PropAir(p_out,T_out)
    Cratio_air=Cp_air/Cv_air
    Area=0.002
    p_amb=101000
    if switch==1 
      k_out=0.3629*10^-5
      flow_out=k_out*(p_out-p_next)
    #   elseif (p_amb/p_out)>((2/(Cratio_air+1))^(1-1/Cratio_air))
    #       flow_out=Cp_air*Area*p_out*(p_amb/p_out)^(1/Cratio_air)*sqrt(2*Cratio_air*(1-(p_amb/p_out)^(1-1/Cratio_air))/(Cratio_air-1))/(sqrt(R_air*T_out))
    #         else
    #       flow_out=Cp_air*Area*p_out*sqrt(Cratio_air)*(2/(Cratio_air+1)^((Cratio_air+1)/(2*(Cratio_air-1))))/(sqrt(R_air*T_out))
    #   end
    end
    return(flow_out,rho_air,R_air)
end

(flow_out,rho_air,R_air)=pipeout(1,T_out,p_out,p_next)
volu=0.02
eqs=(D(m)~flow_in-flow_out,D(p_out)~Cratio*R_air*(flow_in*T_in-flow_out*T_out)/volu)
sm_air=ODESystem(eqs,name=:sm_air)
v_sm=[m=>rho*volu,p_out=>600000]
p_sm=[flow_in=>1.5*10^-4,T_out=>355.15,p_next=>600000]
prob_sm=ODEProblem(sm_air,v_sm,(0,100),p_sm)
sol_sm=solve(prob_sm,Rodas5())

(flow_out,rho)=pipeout(2,T_out,p_out,p_next)
volu=0.005
eqs=(D(m)~flow_in-flow_out,D(p_out)~Cratio*R_air*(flow_in*T_in-flow_out*T_out)/volu)
rm_air=ODESystem(eqs,name=:rm_air)
v_rm=[m=>rho*volu_rm,p_out=>600000]
p_rm=[flow_in=>1.5*10^-4,T_out=>355.15,p_next=>101000]
prob_rm=ODEProblem(rm_air,v_rm,(0,100),p_rm)
sol_rm=solve(prob_rm,Rodas5())