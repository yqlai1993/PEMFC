using ModelingToolkit, DifferentialEquations
using PyCall
using Plots
@pyimport CoolProp.CoolProp as CP
@variables t

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
    elseif medium == 1 && p==0 #"WL"
        R=8.3145/0.018
        rho=-3.9837*10^-3*T^2+2.1274*T+7.1666*10^2
        cv=-3.0788*10^-4*T^3+0.30602*T^2-1.0086*10^2*T+1.5207*10^4
        cp=cv
        k=-7.0674*10^-6*T^2+5.6853*10^-3*T-0.45602
        mu=-3.7024*10^-9*T^3+3.7182*10^-6*T^2-1.2513*10^-3*T+0.14156
    elseif medium == 0 #"steel"
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
    return(p)
end
@register psat(T)
function latent(T) "汽化潜热"
    T=T-273.15
    H=(45070-41.9*T+3.44*10^-3*T^2+2.54*10^-6*T^3-8.98*10^-10*T^4)/0.018
    return(H)
end
@register latent(T)
function stateEq(m,R,T,volu) "状态方程"
    p=m*R*T/volu
    return(p)
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
        println("mvapliq error")
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
    if cv == 0 #阴极
        M1=M_O2
        M2=M_N2
    elseif cv == 1 #阳极
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
        println("memW error")
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
function hcoeff(width,height,length,p,Ts,Ta,medium,condition)
    area=width*height
    D=2*area/(width+height)
    T_ave=(Ts+Ta)/2
    vel=(flow_in+flow_out)/area/2
    ρ=Prop(p,T_ave,medium,1)
    μ=Prop(p,T_ave,medium,5)
    ν=μ/ρ
    Re=ρ*vel*D/μ
    λ=Prop(p,T_ave,medium,4)
    Pr=μ*cp/λ
    Gr=9.81*(1/Ts)*(Ts-Ta)*length^3/ν^2
    Ra=Gr*Pr
    if condition==0 #"ExForConv"
        C=0.102
        n=0.675
        m=1/3
        Nu=C*Re^n*Pr^m
    elseif condition==1 #"Hor&NatConv"
        if Ra<10^7
            C=0.54
        else
            C=0.15
        end
        n=1/4
        Nu=C*Ra^n
    elseif condition==2 #"Ver&NatConv"
        Nu=(0.825+0.387*Ra^(1/6)/(1+(0.5/Pr)^(9/16))^(8/27))^2
    else
        println("hcoeff error")
    end
    h=Nu*λ/D 
    return(h) 
end
@register hcoeff(width,height,length,p,Ts,Ta,medium,condition)
function fcoeff(width,height,length,Num_bend,p,Ts,Ta,medium)
    area=width*height
    D=2*area/(width+height)
    T_ave=(Ts+Ta)/2
    vel=(flow_in+flow_out)/area/2
    ρ=Prop(p,T_ave,medium,1)
    μ=Prop(p,T_ave,medium,5)
    Re=ρ*vel*D/μ
    if Re<2000
        f=(55+41.5exp(-3.4/(width/height)))/Re
    else
        f=1/(1.14-2*log(0.1))
    end
    k=30*f
    dp=(f*length/D+Num_bend*k)*ρ*vel^2/2
    return(dp)
end
@register fcoeff(width,height,length,Num_bend,p,Ts,Ta,medium)
function ini(T,p,volu,c,output)
    m=Prop(T,p,c,1)*volu
    if output==1
        return(T)
    elseif output==2
        return(p)
    elseif output==3
        return(m)
    else
        println("ini error")
    end
end
@register ini(T,p,volu,c,output)

function passage(;name,volu,mix,g1,g2,cv,y_g1_in,h,m_cell) "流道"
    vs= @variables m_g1(t)=101325*0.01/259.9/298.15/120 m_g2(t)=0 m_water(t)=0 T(t)=298.15 T_in(t) p_water_in(t) p_mix_in(t) p_in(t) Rhum_in(t) m_wv(t) m_wl(t) p_g1(t) p_g2(t) p_water(t) p_mix(t) p(t) Rhum(t) flow_mix_in(t) flow_in(t)=0.36*10^-3*101325*2 flow_water_in(t) flow_g1_in(t) flow_g2_in(t) flow_out(t)=0.22*10^-4*101325*1 flow_mix_out(t) flow_water_out(t) flow_g1_out(t) flow_g2_out(t) flow_g1_react(t) flow_water_gen(t) flow_water_mem(t) T_cell(t)
    ps= @parameters volu=volu mix=mix g1=g1 g2=g2 cv=cv y_g1_in=y_g1_in h=h m_cell=m_cell
    D=Differential(t)
    eqs=[#入口压力
        p_water_in~Rhum_in*psat(T_in),
        p_mix_in~p_in-p_water_in,
        #内部压力
        m_wv~mvapliq(m_water,psat(T),volu,Prop(1,298.15,1,8),T,1),
        m_wl~mvapliq(m_water,psat(T),volu,Prop(0,298.15,1,8),T,2),
        p_g1~stateEq(m_g1,Prop(101325,298.15,g1,8),T,volu),
        p_g2~stateEq(m_g2,Prop(101325,298.15,g2,8),T,volu),
        p_water~min(stateEq(m_wv,Prop(1,298.15,1,8),T,volu),psat(T)),
        p_mix~p_g1+p_g2,
        p~p_mix+p_water,
        Rhum~p_water/psat(T),
        #入口流量
        flow_mix_in~flowMed(y_g1_in,p_water_in,p_mix_in,flow_in,cv,1),
        flow_water_in~flowMed(y_g1_in,p_water_in,p_mix_in,flow_in,cv,2),
        flow_g1_in~flowMed(y_g1_in,p_water_in,p_mix_in,flow_in,cv,3),
        flow_g2_in~flowMed(y_g1_in,p_water_in,p_mix_in,flow_in,cv,4),
        #出口流量
        flow_mix_out~flowMed(p_g1/p_mix,p_water,p_mix,flow_out,cv,1),
        flow_water_out~flowMed(p_g1/p_mix,p_water,p_mix,flow_out,cv,2),
        flow_g1_out~flowMed(p_g1/p_mix,p_water,p_mix,flow_out,cv,3),
        flow_g2_out~flowMed(p_g1/p_mix,p_water,p_mix,flow_out,cv,4),
        #微分方程
        D(m_g1)~flow_g1_in-flow_g1_out-flow_g1_react,
        D(m_g2)~flow_g2_in-flow_g2_out,
        D(m_water)~flow_water_in-flow_water_out+flow_water_gen+flow_water_mem,
        D(T)~(h*(T_cell-T)+(Prop(p_mix_in,T_in,mix,3)*flow_mix_in)*T_in-(Prop(p_g1,T,g1,3)*(flow_g1_out+flow_g1_react)+Prop(p_g2,T,g2,3)*flow_g2_out)*T)/(Prop(p_g1,T,g1,3)*m_g1+Prop(p_g2,T,g2,3)*m_g2)
        ]
    ODESystem(eqs,t,vs,ps;name=name)
end
# @named model=passage(volu=0.001,mix=2,g1=3,g2=4,cv=0,h=2,m_stack=10)
# model_simple=structural_simplify(model)
# var0=[m_g1=>101325*0.01/259.28/298.15,m_g2=>0,m_water=>0,T_in=>298.15,T_cell=>298.15,p_in=>101325*3,flow_in=>0.36*10^-4*2*101325,flow_out=>0.22*10^-4*101325,flow_g1_react=>0,flow_water_gen=>0,flow_water_mem=>0]
# tspan=(0,10)
# prob=ODEProblem(model_simple,var0,tspan)
# sol=solve(prob,Rodas5())
# println(sol)

function membrane(;name,area,thick,m_cell,h_ca,h_an,h_amb,T_amb,dH_react) "质子交换膜"
    vs= @variables T(t)=298.15 current(t)=10 flow_H2_react(t) flow_O2_react(t) flow_water_gen(t) flow_water_mem(t) Rhum_ca(t) Rhum_an(t) volt(t) p_H2(t) p_O2(t) T_ca(t) T_an(t)
    ps= @parameters  area=area thick=thick m_cell=m_cell h_ca=h_ca h_an=h_an h_amb=h_amb T_amb=T_amb dH_react=dH_react
    D=Differential(t)
    eqs=[
        flow_H2_react~flowRct(current,1),
        flow_O2_react~flowRct(current,2),
        flow_water_gen~flowRct(current,3),
        flow_water_mem~flowMem(Rhum_ca,Rhum_an,T,current,area,thick),
        volt~voltage(T,p_H2,p_O2,current),
        D(T)~(h_ca*(T_ca-T)+h_an*(T_an-T)+h_amb*(T_amb-T)+flow_H2_react*dH_react-volt*current)/(Prop(101325,298.15,0,3)*m_cell),
        D(current)~0
        ]
    ODESystem(eqs,t,vs,ps;name=name)
end

function cell(;name) "单电池"
    @named ca=passage(volu=0.010/120,mix=2,g1=3,g2=4,cv=0,y_g1_in=0.21,h=2/120,m_cell=10/120)
    @named an=passage(volu=0.005/120,mix=5,g1=5,g2=4,cv=1,y_g1_in=1.00,h=2/120,m_cell=10/120)
    @named mem=membrane(area=232*10^-4,thick=1.275*10^-4,m_cell=10/120,h_ca=2/120,h_an=2/120,h_amb=17/120,T_amb=298.15,dH_react=1.196E+8)
    vs= @variables m_g1_c(t) m_g1_a(t) m_g2_c(t) m_g2_a(t) m_water_c(t) m_water_a(t) p_ci(t) p_ai(t) p_c(t) p_a(t) T_ci(t) T_ai(t) T_c(t) T_a(t) T(t) flow_ci(t) flow_ai(t) flow_co(t) flow_ao(t) Rhum_ci(t) Rhum_ai(t) Rhum_c(t) Rhum_a(t) current(t) volt(t)
    eqs=[
        m_g1_c~ca.m_g1,
        m_g1_a~an.m_g1,
        m_g2_c~ca.m_g2,
        m_g2_a~an.m_g2,
        m_water_c~ca.m_water,
        m_water_a~an.m_water,
        p_ci~ca.p_in,
        p_ai~an.p_in,
        p_c~ca.p,
        p_a~an.p,
        T_ci~ca.T_in,
        T_ai~an.T_in,
        T_c~ca.T,
        T_a~an.T,
        T~mem.T,
        flow_ci~ca.flow_in,
        flow_ai~an.flow_in,
        flow_co~ca.flow_out,
        flow_ao~an.flow_out,
        Rhum_ci~ca.Rhum_in,
        Rhum_ai~an.Rhum_in,
        Rhum_c~ca.Rhum,
        Rhum_a~an.Rhum,
        current~mem.current,
        volt~mem.volt,
        0~ca.p_g1-mem.p_O2,
        0~an.p_g1-mem.p_H2,
        0~ca.T-mem.T_ca,
        0~an.T-mem.T_an,
        0~ca.T_cell-mem.T,
        0~an.T_cell-mem.T,
        0~ca.flow_g1_react-mem.flow_O2_react,
        0~an.flow_g1_react-mem.flow_H2_react,
        0~ca.flow_water_gen-mem.flow_water_gen,
        0~an.flow_water_gen,
        0~ca.flow_water_mem-mem.flow_water_mem,
        0~an.flow_water_mem+mem.flow_water_mem,
        0~ca.Rhum-mem.Rhum_ca,
        0~an.Rhum-mem.Rhum_an
        ]
    compose(ODESystem(eqs,t,vs,[];name=name),ca,an,mem)
end
# @named model=cell()
# model_simple=structural_simplify(model)
# var0=[p_ci=>101325*3,p_ai=>101325*3,p_co=>101325,p_ao=>101325,T_ci=>298.15,T_ai=>298.15,T_co=>298.15,T_ao=>298.15,model.flow_ci=>101325*2*0.36*10^-4,model.flow_ai=>101325*2*0.36*10^-4,model.flow_co=>101325*0.22*10^-4,model.flow_ao=>101325*0.22*10^-4]
# tspan=(0,10)
# prob=ODEProblem(model_simple,var0,tspan)
# sol=solve(prob,Rodas5())
# println(sol)

function stack(;name,num) "电池堆"
    @named c1=cell()
    
    vs= @variables m_g1_c(t) m_g1_a(t) m_g2_c(t) m_g2_a(t) m_water_c(t) m_water_a(t) p_ci(t) p_ai(t) p_c(t) p_a(t) T_ci(t) T_ai(t) T_c(t) T_a(t) T_c1(t) Rhum_ci(t) Rhum_ai(t) Rhum_c(t) Rhum_a(t) current(t) volt(t) flow_ci(t) flow_ai(t) flow_co(t) flow_ao(t)
    ps= @parameters num=num
    eqs=[
        m_g1_c~c1.m_g1_c*num,
        m_g1_a~c1.m_g1_a*num,
        m_g2_c~c1.m_g2_c*num,
        m_g2_a~c1.m_g2_a*num,
        m_water_c~c1.m_water_c*num,
        m_water_a~c1.m_water_a*num,
        p_ci~c1.p_ci,
        p_ai~c1.p_ai,
        p_c~c1.p_c,
        p_a~c1.p_a,
        T_ci~c1.T_ci,
        T_ai~c1.T_ai,
        T_c~c1.T_c,
        T_a~c1.T_a,
        T_c1~c1.T,
        flow_ci~c1.flow_ci*num,
        flow_ai~c1.flow_ai*num,
        flow_co~c1.flow_co*num,
        flow_ao~c1.flow_ao*num,
        Rhum_ci~c1.Rhum_ci,
        Rhum_ai~c1.Rhum_ai,
        Rhum_c~c1.Rhum_c,
        Rhum_a~c1.Rhum_a,
        current~c1.current,
        volt~c1.volt*num
        ]
    compose(ODESystem(eqs,t,vs,ps;name=name),c1)
end


function up(;name) "上游部件"
    vs= @variables p(t)=101325 T(t)=298.15 Rhum(t)=0.5
    D=Differential(t)
    eqs=[D(p)~0,D(T)~0,D(Rhum)~0]
    ODESystem(eqs,t,vs,[];name=name)
end

function down(;name) "下游部件"
    vs= @variables p(t)=101325 T(t)
    D=Differential(t)
    eqs=[D(p)~0]
    ODESystem(eqs,t,vs,[];name=name)
end

@named uc=up()
@named ua=up()
@named dc=down()
@named da=down()
@named st=stack(num=120)
@variables t
@parameters k_ca_in k_ca_out k_an_in k_an_out
eqs=[
    0~st.T_ci-uc.T,
    0~st.T_ai-ua.T,
    0~st.T_c-dc.T,
    0~st.T_a-da.T,
    0~st.p_ci-uc.p,
    0~st.p_ai-ua.p,
    0~st.flow_ci-flowtot(k_ca_in,uc.p,st.p_c),
    0~st.flow_ai-flowtot(k_an_in,ua.p,st.p_a),
    0~st.flow_co-flowtot(k_ca_out,dc.p,st.p_c),
    0~st.flow_ao-flowtot(k_an_out,da.p,st.p_a),
    0~st.Rhum_ci-uc.Rhum,
    0~st.Rhum_ai-ua.Rhum
    ]
@named fc=ODESystem(eqs,t)
@named model=compose(fc,[uc,ua,dc,da,st])
model_simple=structural_simplify(model)

var0=[uc.p=>101325*5,ua.p=>101325*7,dc.p=>101325,da.p=>101325,st.c1.an.m_g1=>101325*0.005/4157.25/298.15/120]
tspan=(0.0,10.0)
par0=[k_ca_in=>0.36*10^-3,k_ca_out=>0.22*10^-5,k_an_in=>0.36*10^-5,k_an_out=>0.22*10^-5]
prob=ODEProblem(model_simple,var0,tspan,par0)
sol=solve(prob,Rodas4(),saveat=1,maxiters=1e7)

using DataFrames
df=DataFrame(sol)
println(df)