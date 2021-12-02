using ModelingToolkit,OrdinaryDiffEq
# using Unitful,RecursiveArrayTools

function Prop(p,T,medium,output)
    if medium == 2 #"Air"
        R=286.987
        rho=p/(R*T)
        cv=R*(3.653-1.337e-3*T+3.294e-6*T^2-1.913e-9*T^3+0.2763e-12*T^4-1)
        cp=cv+R
        k=0.009748221+5.27354*10^-10*p+5.5243*10^-5*T
        mu=6.93093*10^-6+2.35465*10^-13*p+3.85177*10^-8*T
    elseif medium == 5 #"H2"
        R=8.3145/0.002 
        rho=p/(R*T)
        cp=(6.88+0.000066*T+0.279*10^-6*T^2)/2*4.1868*1000
        cv=cp-R
    elseif medium == 4 #"N2"
        R=8.3145/0.028
        rho=p/(R*T)
        cp=(6.7619+0.000871*T-1.643*10^-7*T^2)/28*4.1868*1000
        #cp=(6.3+0.001819*T-0.345*10^-6*T^2)/28*4.1868*1000
        cv=cp-R
    elseif medium == 3 #"O2"
        R=8.3145/0.032
        rho=p/(R*T)
        cp=(6.9647+1.0971*10^-3*T-2.03*10^-7*T^2)/32*4.1868*1000
        #cp=(7.025+0.00185*T-abs((T-300)*0.075/300))/32*4.1868*1000
        cv=cp-R
    elseif medium == 1 && p==1 #"WV"
            R=8.3145/0.018 
            rho=p/(R*T)
            cp=(7.75+0.0016*T-1.6*10^-7*T^2)/18*4.1868*1000
            cv=cp-R 
    elseif medium == 1 && p == 0  #"WL"
            R=8.3145/0.018
            rho=-3.9837*10^-3*T^2+2.1274*T+7.1666*10^2
            cv=-3.0788*10^-4*T^3+0.30602*T^2-1.0086*10^2*T+1.5207*10^4
            cp=cv
            k=-7.0674*10^-6*T^2+5.6853*10^-3*T-0.45602
            mu=-3.7024*10^-9*T^3+3.7182*10^-6*T^2-1.2513*10^-3*T+0.14156
     
    elseif medium == 0 #"steel"
        R=8.3145/0.056
        cp=350     
    else
        println("Prop medium error")
    end
    
    if output==1
        return(rho)
    elseif output==2
        return(cv)
    elseif output==3
        return(cp)
    elseif output==8
        return(R)
    else
        println("Prop output error")
    end
end         
@register Prop(p,T,medium,output)

function psat(T) "饱和水蒸汽压"
    T=T-273.15
    # a=-2.1794+0.02953*T-9.1837*10^-5*T^2+1.4454*10^-7*T^3
    a=-2.18+0.0295*T-9.18*10^-5*T^2+1.44*10^-7*T^3
    p=10^(a+5)
end
@register psat(T)

function latent(T) "汽化潜热"
    T=T-273.15
    H=(45070-41.9*T+3.44*10^-3*T^2+2.54*10^-6*T^3-8.98*10^-10*T^4)/0.018
end
@register latent(T)

function stateEq(m,R,T,volu) "状态方程"
    p=m*R*T/volu
end
@register stateEq(m,R,T,volu)

function mvapliq(m_water,p_sat,volu,R_water,T,output) "气态和液态质量"
    m_water_max=p_sat*volu/R_water/T
    if m_water<=m_water_max
       m_water_vap=m_water
       m_water_liq=0
     else
       m_water_vap=m_water_max
       m_water_liq=m_water-m_water_max 
    end
    if output==1
        return(m_water_vap)
    elseif output==2
        return(m_water_liq)
    else
        println("mvapliq output error")
    end
end
@register mvapliq(m_water,p_sat,volu,R_water,T,output)

function flowtot(k,p_in,p_out) "进出口总流量"
    flow=k*(p_in-p_out)
    return(flow)
end
@register flowtot(k,p_in,p_out)

function flowMed(y1,p_water,p_gas,flow_tot,cv,output) "工质质量流量"
    M_H2=2*10^-3
    M_CO=28*10^-3
    M_O2=32*10^-3
    M_N2=28*10^-3
    M_water=18*10^-3
    if cv == 0 
        M1=M_O2
        M2=M_N2
    elseif cv == 1 
        M1=M_H2
        M2=M_N2
    else
        println("flowMed cv error")
    end  
    y2=1-y1
    M_gas=y1*M1+y2*M2
    omega=M_water*p_water/M_gas/p_gas
    flow_gas=flow_tot/(1+omega)
    flow_water=flow_tot-flow_gas
    flow_g1=flow_gas*y1*M1/M_gas
    flow_g2=flow_gas*y2*M2/M_gas
    if output==1
        return(flow_gas)
    elseif output==2
        return(flow_water)
    elseif output==3
        return(flow_g1)
    elseif output==4
        return(flow_g2)
    else
        println("flowMed output error")
    end
end
@register flowMed(y1,p_water,p_gas,flow_tot,cv,output)

function memW(a,output) "膜含水量"
    Den_mem=2000
    M_mem=1.1
    if a<=1
        content_w=0.043+17.81*a-39.85*a^2+36*a^3
    elseif a<=3
        content_w=14+1.4*(a-1)
    else
        content_w=16.8
    end
    MolarDen_w=Den_mem*content_w/M_mem

    if output==1
        return(content_w)
    elseif output==2
        return(MolarDen_w)
    else
        println("memW output error")
    end
end
@register memW(a,output)

function flowMem(Rhum_ca,Rhum_an,T_cell,current,area_mem,thick_mem) "膜水合流量"
    M_water=0.018
    F=96485
    σ_ca=memW(Rhum_ca,1)
    c_ca=memW(Rhum_ca,2)
    σ_an=memW(Rhum_an,1)
    c_an=memW(Rhum_an,2)
    σ_ave=(σ_an+σ_ca)/2
    n_d=0.0029*(σ_ca+σ_an)^2/4+0.05*σ_ave-3.4e-19
    D_w=1.25e-6*exp(2416*(1/303-1/T_cell))
    currentDen=current/area_mem
    N_water_mem=n_d*currentDen/F-D_w*(c_ca-c_an)/thick_mem
    flow_water_mem=N_water_mem*M_water*area_mem
    return(flow_water_mem)
end
@register flowMem(Rhum_ca,Rhum_an,T_cell,current,area_mem,thick_mem)

function flowRct(current,output) "反应过程量"
    F=96485
    M_H2=2*10^-3
    M_O2=32*10^-3
    M_water=18*10^-3
    flow_H2_react=M_H2*current/2/F
    flow_O2_react=M_O2*current/4/F
    flow_water_gen=M_water*current/2/F
    if output==1
        return(flow_H2_react)
    elseif output==2
        return(flow_O2_react)
    elseif output==3
        return(flow_water_gen)
    else
        println("flowRct error")
    end
end
@register flowRct(current,output)

function voltage(T_cell,p_H2_an,p_O2_ca,current) "电压"
    R=8.3145
    F=96485
    volt_nst=1.229-8.5e-4*(T_cell-298.15)+R*T_cell*(log(p_H2_an)+0.5*log(p_O2_ca))/(2*F)
    volt_act=0.9514-3.12e-3*T_cell+1.87e-4*T_cell*log(current)-7.4e-5*T_cell*log(p_O2_ca/(5.08e6*exp(-498/T_cell)))
    #resistivity_mem=(thick_mem/area_mem)*181.6*(1+0.03*currentDen+0.062*(T_cell/303)^2*currentDen^2.5)/((sigma_an-0.634-3*currentDen)*exp(4.18*(1-303/T_cell)))
    resistivity_mem=0.01605-3.5e-5*T_cell+8e-5*current
    volt_ohm=current*resistivity_mem
    volt_cell=volt_nst-volt_act-volt_ohm
    return(volt_cell)
end 
@register voltage(T_cell,p_H2_an,p_O2_ca,current)

function iniPara(output)
    T_in=298.15
    p_in=101325
    Rhum_in=0.5
    y_ca_in=0.21
    y_an_in=1
    if output==1
        return(T_in)
    elseif output==2
        return(p_in)
    elseif output==3
        return(Rhum_in)
    elseif output==4
        return(y_ca_in)
    elseif output==5
        return(y_an_in)
    else
        println("iniPara error")
    end
end
@register iniPara(output)

function constPara(output)
    T_amb=298.15
    p_amb=101325
    thick_mem=1.275*10^-4
    area_mem=232*10^-4
    m_stack=10
    num_cell=120
    volu_ca=0.01
    volu_an=0.005
    k_in=0.36*10^-5
    k_out=0.22*10^-5
    if output==1
        return(T_amb)
    elseif output==2
        return(p_amb)
    elseif output==3
        return(thick_mem)
    elseif output==4
        return(area_mem)
    elseif output==5
        return(m_stack)
    elseif output==6
        return(num_cell)
    elseif output==7
        return(volu_ca)
    elseif output==8
        return(volu_an)
    elseif output==9
        return(k_in)
    elseif output==10
        return(k_out)
    else
        println("constPara error")
    end
end
@register constPara(output)

@parameters dH_react,current,h_amb,h_ca,h_an
@variables t,m_O2_ca(t),m_N2_ca(t),m_water_ca(t),m_H2_an(t),m_N2_an(t),m_water_an(t),T_cell(t),T_ca(t),T_an(t),
    p_water_ca_in(t),p_air_ca_in(t),Cp_O2_ca_in(t),Cp_N2_ca_in(t),Cp_wv_ca_in(t),Cv_wl_ca_in(t),
    p_sat_ca(t),m_wv_ca(t),m_wl_ca(t),p_O2_ca(t),p_N2_ca(t),p_water_ca(t),p_air_ca(t),p_ca(t),Rhum_ca(t),Cp_O2_ca(t),Cp_N2_ca(t),Cp_wv_ca(t),Cv_wl_ca(t),
    flow_ca_in(t),flow_air_ca_in(t),flow_water_ca_in(t),flow_O2_ca_in(t),flow_N2_ca_in(t),
    flow_ca_out(t),flow_air_ca_out(t),flow_water_ca_out(t),flow_O2_ca_out(t),flow_N2_ca_out(t),
    p_water_an_in(t),p_mix_an_in(t),Cp_H2_an_in(t),Cp_N2_an_in(t),Cp_wv_an_in(t),Cv_wl_an_in(t),
    p_sat_an(t),m_wv_an(t),m_wl_an(t),p_H2_an(t),p_N2_an(t),p_water_an(t),p_mix_an(t),p_an(t),Rhum_an(t),Cp_H2_an(t),Cp_N2_an(t),Cp_wv_an(t),Cv_wl_an(t),
    flow_an_in(t),flow_mix_an_in(t),flow_water_an_in(t),flow_H2_an_in(t),flow_N2_an_in(t),
    flow_an_out(t),flow_mix_an_out(t),flow_water_an_out(t),flow_H2_an_out(t),flow_N2_an_out(t),
    flow_O2_react(t),flow_H2_react(t),flow_water_gen(t),
    σ_ca(t),c_ca(t),σ_an(t),c_an(t),σ_ave(t),n_d(t),D_w(t),currentDen(t),N_water_mem(t),flow_water_mem(t),
    volt_cell(t),power_cell(t),eff_cell(t)
D=Differential(t)

eqs=[#阴极入口压力
    p_water_ca_in~iniPara(3)*psat(iniPara(1)),
    p_air_ca_in~iniPara(2)-p_water_ca_in,    
    Cp_wv_ca_in~Prop(1,iniPara(1),1,3),
    Cv_wl_ca_in~Prop(0,iniPara(1),1,2),
    #阴极内部压力
    p_sat_ca~psat(T_ca),
    m_wv_ca~mvapliq(m_water_ca,p_sat_ca,constPara(7),Prop(1,298.15,1,8),T_ca,1),
    m_wl_ca~mvapliq(m_water_ca,p_sat_ca,constPara(7),Prop(0,298.15,1,8),T_ca,2),
    p_O2_ca~stateEq(m_O2_ca,Prop(101325,298.15,3,8),T_ca,constPara(7)),
    p_N2_ca~stateEq(m_N2_ca,Prop(101325,298.15,4,8),T_ca,constPara(7)),
    p_water_ca~stateEq(m_wv_ca,Prop(1,298.15,1,8),T_ca,constPara(7)), 
    p_air_ca~p_O2_ca+p_N2_ca,
    p_ca~p_air_ca+p_water_ca,
    Rhum_ca~p_water_ca/p_sat_ca,
    Cp_O2_ca~Prop(p_O2_ca,T_ca,3,3),
    Cp_N2_ca~Prop(p_N2_ca,T_ca,4,3),
    Cp_wv_ca~Prop(1,T_ca,1,3),
    Cv_wl_ca~Prop(0,T_ca,1,2),
    #阴极入口流量
    flow_ca_in~flowtot(constPara(9)*100,iniPara(2)*5,p_ca),
    flow_air_ca_in~flowMed(iniPara(4),p_water_ca_in,p_air_ca_in,flow_ca_in,0,1),
    flow_water_ca_in~flowMed(iniPara(4),p_water_ca_in,p_air_ca_in,flow_ca_in,0,2),
    flow_O2_ca_in~flowMed(iniPara(4),p_water_ca_in,p_air_ca_in,flow_ca_in,0,3),
    flow_N2_ca_in~flowMed(iniPara(4),p_water_ca_in,p_air_ca_in,flow_ca_in,0,4),
    Cp_O2_ca_in~Prop(p_air_ca_in*flow_O2_ca_in/0.032/(flow_O2_ca_in/0.032+flow_N2_ca_in/0.028),298.15,3,3),
    Cp_N2_ca_in~Prop(p_air_ca_in*flow_N2_ca_in/0.028/(flow_O2_ca_in/0.032+flow_N2_ca_in/0.028),298.15,4,3),
    #阴极出口流量
    flow_ca_out~flowtot(constPara(10),p_ca,constPara(2)),
    flow_air_ca_out~flowMed(p_O2_ca/p_air_ca,p_water_ca,p_air_ca,flow_ca_out,0,1),
    flow_water_ca_out~flowMed(p_O2_ca/p_air_ca,p_water_ca,p_air_ca,flow_ca_out,0,2),
    flow_O2_ca_out~flowMed(p_O2_ca/p_air_ca,p_water_ca,p_air_ca,flow_ca_out,0,3),
    flow_N2_ca_out~flowMed(p_O2_ca/p_air_ca,p_water_ca,p_air_ca,flow_ca_out,0,4),
    #阳极入口压力
    p_water_an_in~iniPara(3)*psat(iniPara(1)),
    p_mix_an_in~iniPara(2)-p_water_an_in,
    Cp_wv_an_in~Prop(1,iniPara(1),1,3),
    Cv_wl_an_in~Prop(0,iniPara(1),1,2),
    #阳极内部压力
    p_sat_an~psat(T_an),
    m_wv_an~mvapliq(m_water_an,p_sat_an,constPara(8),Prop(1,298.15,1,8),T_an,1),
    m_wl_an~mvapliq(m_water_an,p_sat_an,constPara(8),Prop(0,298.15,1,8),T_an,2),
    p_H2_an~stateEq(m_H2_an,Prop(101325,298.15,5,8),T_an,constPara(8)),
    p_N2_an~stateEq(m_N2_an,Prop(101325,298.15,4,8),T_an,constPara(8)),
    p_water_an~stateEq(m_wv_an,Prop(1,298.15,1,8),T_an,constPara(8)),
    p_mix_an~p_H2_an+p_N2_an,
    p_an~p_mix_an+p_water_an,
    Rhum_an~p_water_an/p_sat_an,
    Cp_H2_an~Prop(p_H2_an,T_an,5,3),
    Cp_N2_an~Prop(p_N2_an,T_an,4,3),
    Cp_wv_an~Prop(1,T_an,1,3),
    Cv_wl_an~Prop(0,T_an,1,2),
    #阳极入口流量
    flow_an_in~flowtot(constPara(9),iniPara(2)*7,p_an),
    flow_mix_an_in~flowMed(iniPara(5),p_water_an_in,p_mix_an_in,flow_an_in,1,1),
    flow_water_an_in~flowMed(iniPara(5),p_water_an_in,p_mix_an_in,flow_an_in,1,2),
    flow_H2_an_in~flowMed(iniPara(5),p_water_an_in,p_mix_an_in,flow_an_in,1,3),
    flow_N2_an_in~flowMed(iniPara(5),p_water_an_in,p_mix_an_in,flow_an_in,1,4),
    Cp_H2_an_in~Prop(p_mix_an_in*flow_H2_an_in/0.002/(flow_H2_an_in/0.002+flow_N2_an_in/0.028),298.15,5,3),
    Cp_N2_an_in~Prop(p_mix_an_in*flow_N2_an_in/0.028/(flow_H2_an_in/0.002+flow_N2_an_in/0.028),298.15,4,3),
    #阳极出口流量
    flow_an_out~flowtot(constPara(10),p_an,constPara(2)),
    flow_mix_an_out~flowMed(1,p_water_an,p_mix_an,flow_an_out,1,1),
    flow_water_an_out~flowMed(1,p_water_an,p_mix_an,flow_an_out,1,2),
    flow_H2_an_out~flowMed(1,p_water_an,p_mix_an,flow_an_out,1,3),
    flow_N2_an_out~flowMed(1,p_water_an,p_mix_an,flow_an_out,1,4),
    #反应量和生成量
    flow_H2_react~constPara(6)*flowRct(current,1),
    flow_O2_react~constPara(6)*flowRct(current,2),
    flow_water_gen~constPara(6)*flowRct(current,3),
    #膜水合流量
    flow_water_mem~constPara(6)*flowMem(Rhum_ca,Rhum_an,T_cell,current,constPara(4),constPara(3)),
    #电压
    volt_cell~voltage(T_cell,p_H2_an,p_O2_ca,current),
    power_cell~volt_cell*current,
    T_ca~T_cell,
    T_an~T_cell,
    #控制方程
    D(m_O2_ca)~flow_O2_ca_in-flow_O2_ca_out-flow_O2_react,
    D(m_N2_ca)~flow_N2_ca_in-flow_N2_ca_out,
    D(m_water_ca)~flow_water_ca_in-flow_water_ca_out+flow_water_gen+flow_water_mem,
    D(m_H2_an)~flow_H2_an_in-flow_H2_an_out-flow_H2_react,
    D(m_N2_an)~flow_N2_an_in-flow_N2_an_out,
    D(m_water_an)~flow_water_an_in-flow_water_an_out-flow_water_mem,
    # D(T_ca)~(h_ca*(T_cell-T_ca)+(Cp_O2_ca_in*flow_O2_ca_in+Cp_N2_ca_in*flow_N2_ca_in+Cp_wv_ca_in*flow_water_ca_in)*iniPara(1)-(Cp_O2_ca*(flow_O2_ca_out+flow_O2_react)+Cp_N2_ca*flow_N2_ca_out+Cp_wv_ca*flow_water_ca_out)*T_ca)/(Cp_O2_ca*m_O2_ca+Cp_N2_ca*m_N2_ca+Cp_wv_ca*m_wv_ca+Cv_wl_ca*m_wl_ca),
    # D(T_an)~(h_an*(T_cell-T_an)+(Cp_H2_an_in*flow_H2_an_in+Cp_N2_an_in*flow_N2_an_in+Cp_wv_an_in*flow_water_an_in)*iniPara(1)-(Cp_H2_an*(flow_H2_an_out+flow_H2_react)+Cp_N2_an*flow_N2_an_out+Cp_wv_ca*flow_water_ca_out)*T_an)/(Cp_H2_an*m_H2_an+Cp_N2_an*m_N2_an+Cp_wv_an*m_wv_an+Cv_wl_an*m_wl_an),
    D(T_cell)~(h_ca*(T_ca-T_cell)+h_an*(T_an-T_cell)+h_amb*(constPara(1)-T_cell)+flow_H2_react*dH_react-constPara(6)*power_cell
                    +(Cp_O2_ca_in*flow_O2_ca_in+Cp_N2_ca_in*flow_N2_ca_in+Cp_wv_ca_in*flow_water_ca_in)*iniPara(1)-(Cp_O2_ca*(flow_O2_ca_out+flow_O2_react)+Cp_N2_ca*flow_N2_ca_out+Cp_wv_ca*flow_water_ca_out)*T_ca
                    +(Cp_H2_an_in*flow_H2_an_in+Cp_N2_an_in*flow_N2_an_in+Cp_wv_an_in*flow_water_an_in)*iniPara(1)-(Cp_H2_an*(flow_H2_an_out+flow_H2_react)+Cp_N2_an*flow_N2_an_out+Cp_wv_ca*flow_water_ca_out)*T_an)/(Prop(101000,298.15,0,3)*constPara(5)+Cp_O2_ca*m_O2_ca+Cp_N2_ca*m_N2_ca+Cp_wv_ca*m_wv_ca+Cv_wl_ca*m_wl_ca+Cp_H2_an*m_H2_an+Cp_N2_an*m_N2_an+Cp_wv_an*m_wv_an+Cv_wl_an*m_wl_an)
    ]
@named cell0=ODESystem(eqs,t)
cell_simple=structural_simplify(cell0)

var=[m_O2_ca=>101325*0.01/259.28/298.15,m_N2_ca=>0,m_water_ca=>0,m_H2_an=>101325*0.005/4157.25/298.15,m_N2_an=>0,m_water_an=>0,T_ca=>298.15,T_an=>298.15,T_cell=>298.15]
tspan=(0.0,10.0)
par=[dH_react=>1.196e8,current=>10.0,h_amb=>17.0,h_ca=>2.0,h_an=>2.0]
prob=ODEProblem(cell_simple,var,tspan,par)
sol=solve(prob,Rodas4(),saveat=1,maxiters=1e7)

using DataFrames
df=DataFrame(sol)
println(df)

using Plots
p1=plot(sol,vars=[volt_cell])
display(p1)
p2=plot(sol,vars=[m_water_ca,m_water_an])
display(p2)
p3=plot(sol,vars=[m_O2_ca,m_N2_ca,m_H2_an,m_N2_an])
display(p3)
p4=plot(sol,vars=[T_ca,T_an,T_cell])
display(p4)
# plot(p1,p2,p3,p4,layout=(4,1))
p5=plot(sol,vars=[flow_N2_an_in,flow_N2_an_out])
display(p5)
using RecursiveArrayTools
t = collect(range(0,stop=50,length=100)) # 建立时间向量
randomized = VectorOfArray([(sol(t[i])+.5randn(7)) for i in 1:length(t)])
data = convert(Array,randomized)  

using DifferentialEquations,DiffEqParamEstim
cost_function = build_loss_objective(prob,Tsit5(),L2Loss(t,data),maxiters=10000000,verbose=false)
using Optim
result = optimize(cost_function, [1.196e8,20.0,17.0,2.0,4.0])
println(result.minimizer)

