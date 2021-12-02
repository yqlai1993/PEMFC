using Base: sign_mask
println("helloworld")
#for
# s=0
# t=0
# a=[]
# for i in 1:10
#   global t=t-1 
#   global s=s+i
#   push!(a,[s t])
# end
# println(a)



#function
# function aa(x,y)
#   a=2*x+y
#   b=x/3-y/4
#   return c=[a,b]
# end
# (boy,girl)=aa(1,1)
# println(boy)
# println(girl)



#—————————————————————————————————————————— 常微分方程 —————————————————————————————————————————
# using ModelingToolkit
# using OrdinaryDiffEq
# using Plots
# @parameters σ,ρ,β
# @variables t,x(t),y(t),z(t)
# D=Differential(t)
# eqs=[D(x)~σ*(y-x),
#     D(y)~x*(ρ-z)-y,
#     D(z)~x*y-β*z]
# @named sys=ODESystem(eqs,t,[x,y,z],[σ,ρ,β])

# u=[1.0,0.0,0.0]
# tspan=(0.0,6.0)
# p=[10.0,28.0,8/3]
# prob=ODEProblem(sys,u,tspan,p)
# sol=solve(prob,Tsit5())
# println(sol)
# plot(sol,vars=[x,y,z])



#——————————————————————————————————————微分代数方程-中间变量——————————————————————————————————————
# using ModelingToolkit,OrdinaryDiffEq,Plots

# @parameters σ,ρ,β
# @variables t,x(t),y(t),z(t),ω(t),γ(t)
# D=Differential(t)
# eqs=[ω~1/x,
#      γ~ω*2,
#      D(x)~σ*(y-x)*γ,
#      D(y)~x*(ρ-z)-y,
#      D(z)~x*y-β*z]
# @named sys=ODESystem(eqs,t,[x,y,z],[σ,ρ,β])
# sys_simple=structural_simplify(sys)

# u=[x=>1.0,y=>0.0,z=>0.0]
# tspan=(0.0,6.0)
# p=[σ=>10.0,ρ=>28.0,β=>8/3]
# prob=ODEProblem(sys_simple,u,tspan,p)
# sol=solve(prob,Tsit5())
# println(sol)
# println(sol[ω])
# println(sol[γ])
# plot(sol,vars=[x,y,z,ω,γ]) 



#—————————————————————————————————————— 微分代数方程-外部函数——————————————————————————————————————
# using ModelingToolkit,DifferentialEquations,Plots

# @parameters σ,ρ,β
# @variables t x(t) y(t) z(t) ω(t) γ(t) α(t)

# function force(x,switch)
#     if x<=3
#         y=10+x/1000
#         z=sin(x)
#     else
#         y=-10+x/1000
#         z=-sin(x)
#     end
#     if switch==1
#         return(y)
#     else
#         return(z)
#     end
# end
# @register force(x,switch)

# function Force(x,switch) 
#     y=force(x,1)/10
#     z=force(x,0)/0.1
#     if switch==1
#         return(y)
#     else
#         return(z)
#     end
# end
# @register Force(x,switch)

# D=Differential(t)
# eqs=[ω~2*t,
#     γ~1*t,
#     α~Force(t,1),
#     D(x)~ω+σ*(y-x),
#     D(y)~x*(ρ-z)-y+γ,
#     D(z)~x*y-β*z+α]
# @named sys=ODESystem(eqs)
# sys_simple=structural_simplify(sys)

# u=[x=>1.0,y=>0.0,z=>0.0]
# tspan=(0.0,6.0)
# p=[σ=>10.0,ρ=>28.0,β=>8/3]
# prob=ODEProblem(sys_simple,u,tspan,p)
# sol=solve(prob)
# println(sol[γ])
# plot(sol,vars=[x,y,z,ω,γ,α])
  


#———————————————————————————————————————— 基于组件的建模一 ————————————————————————————————————————
# using ModelingToolkit,OrdinaryDiffEq,Plots

# function unit1(;name)
#     @parameters σ,ρ,β
#     @variables t,x(t),y(t),z(t)
#     D=Differential(t)
#     eqs=[D(x)~σ*(y-x),
#         D(y)~x*(ρ-z)-y,
#         D(z)~x*y-β*z]
#     ODESystem(eqs,t,[x,y,z],[σ,ρ,β];name=name)
# end

# function unit2(;name,K)
#     @parameters σ,ρ,β
#     @variables t,x(t),y(t),z(t)
#     D=Differential(t)
#     eqs=[D(x)~σ*(y-x)*K,
#         D(y)~x*(ρ-z)-y,
#         D(z)~x*y-β*z]
#     ODESystem(eqs,t,[x,y,z],[σ,ρ,β];name=name)
# end

# @variables ω(t)
# @named model1=unit1()
# @named model2=unit2(K=2)
# cycle_eqs=[0~model1.x+model2.y+ω]
# @named sys=ODESystem(cycle_eqs,t,systems=[model1,model2])
# sys_simple=structural_simplify(sys)

# u=[model1.x=>1.0,model1.y=>0.0,model1.z=>0.0,
#     model2.x=>2.0,model2.y=>1.0,model2.z=>0.0,
#     ω=>2.0]
# tspan=(0.0,6.0)
# p=[model1.σ=>10.0,model1.ρ=>28.0,model1.β=>8/3,
#     model2.σ=>10.0,model2.ρ=>28.0,model2.β=>8/3]
# prob=ODEProblem(sys_simple,u,tspan,p)
# sol=solve(prob,Tsit5())
# println(sol)
# plot(sol)



#———————————————————————————————————————— 基于组件的建模二 ————————————————————————————————————————
# using ModelingToolkit,OrdinaryDiffEq,Plots
# function v0(output)
#     if output==1
#         x=1
#         return(x)
#     elseif output==2
#         y=0
#         return(y)
#     else
#         println("error")
#     end
# end
# @register v0(a)

# function unit(;name)
#     @parameters σ ρ β
#     vs=@variables x(t)=v0(1) y(t)=v0(2) z(t)=v0(2) α(t)
#     D=Differential(t)
#     eqs=[α~2*x,
#         D(x)~σ*(y-x)+α,
#         D(y)~x*(ρ-z)-y,
#         D(z)~x*y-β*z]
#     # vs=[x,y,z]
#     ps=[σ,ρ,β]     
#     ODESystem(eqs,t,vs,ps;name=name)
# end

# @named c1=unit()
# @named c2=unit()
# @variables t,ω(t)
# cycle_eqs=[0~c1.x+c2.y-ω]
# @named sys=ODESystem(cycle_eqs,t,systems=[c1,c2])
# sys_simple=structural_simplify(sys)

# u0=[c1.x=>1,c1.y=>0,c1.z=>0,c2.x=>2,c2.y=>0,c2.z=>0,ω=>1]
# tspan=(0.0,6.0)
# p0=[c1.σ=>1,c1.ρ=>1,c1.β=>10,c2.σ=>2,c2.ρ=>2,c2.β=>10]
# prob=ODEProblem(sys_simple,u0,tspan,p0)
# sol=solve(prob,Tsit5())
# println(sol)
# plot(sol)



#———————————————————————————————————————— 基于组件的建模三 ————————————————————————————————————————
# using ModelingToolkit, Plots, DifferentialEquations

# @variables t
# @connector function Pin(;name)
#     sts = @variables v(t)=1.0 i(t)=1.0
#     ODESystem(Equation[], t, sts, []; name=name)
# end

# function ModelingToolkit.connect(::Type{Pin}, ps...)
#     eqs = [
#            0 ~ sum(p->p.i, ps) # KCL
#           ]
#     # KVL
#     for i in 1:length(ps)-1
#         push!(eqs, ps[i].v ~ ps[i+1].v)
#     end

#     return eqs
# end

# function Ground(;name)
#     @named g = Pin()
#     eqs = [g.v ~ 0]
#     compose(ODESystem(eqs, t, [], []; name=name), g)
# end

# function OnePort(;name)
#     @named p = Pin()
#     @named n = Pin()
#     sts = @variables v(t)=1.0 i(t)=1.0
#     eqs = [
#            v ~ p.v - n.v
#            0 ~ p.i + n.i
#            i ~ p.i
#           ]
#     compose(ODESystem(eqs, t, sts, []; name=name), p, n)
# end

# function Resistor(;name, R1 = 1.0, R2 = 1.0)
#     @named oneport = OnePort()
#     @unpack v, i = oneport
#     ps = @parameters R1=R1 R2=R2
#     eqs = [
#            v ~ i * (R1+R2)
#           ]
#     extend(ODESystem(eqs, t, [], ps; name=name), oneport)
# end

# function Capacitor(;name, C = 1.0)
#     @named oneport = OnePort()
#     @unpack v, i = oneport
#     ps = @parameters C=C
#     D = Differential(t)
#     eqs = [
#            D(v) ~ i / C
#           ]
#     extend(ODESystem(eqs, t, [], ps; name=name), oneport)
# end

# function ConstantVoltage(;name, V = 1.0)
#     @named oneport = OnePort()
#     @unpack v = oneport
#     ps = @parameters V=V
#     eqs = [
#            V ~ v
#           ]
#     extend(ODESystem(eqs, t, [], ps; name=name), oneport)
# end

# R = 1.0
# C = 1.0
# V = 1.0
# @named resistor = Resistor(R1=R,R2=R)
# @named capacitor = Capacitor(C=C)
# @named source = ConstantVoltage(V=V)
# @named ground = Ground()

# rc_eqs = [
#           connect(source.p, resistor.p)
#           connect(resistor.n, capacitor.p)
#           connect(capacitor.n, source.n, ground.g)
#          ]

# @named _rc_model = ODESystem(rc_eqs, t)
# @named rc_model = compose(_rc_model,[resistor, capacitor, source, ground])
# sys = structural_simplify(rc_model)
# u0 = [
#       capacitor.v => 0.0
#      ]
# prob = ODAEProblem(sys, u0, (0, 10.0))
# sol = solve(prob, Tsit5())
# plot(sol)



#—————————————————————————————————————————— 连续回调 ——————————————————————————————————————————
# using ModelingToolkit
# using DifferentialEquations
# using Plots

# @parameters σ,ρ,β
# @variables t,x(t),y(t),z(t),ω(t)
# D=Differential(t)
# eqs=[D(x)~σ*(y-x)*ω,
#     D(y)~x*(ρ-z)-y,
#     D(z)~x*y-β*z,
#     D(ω)~0]
# @named sys=ODESystem(eqs,t,[x,y,z,ω],[σ,ρ,β])

# u=[x=>1.0,y=>0.0,z=>0.0,ω=>1.0]
# tspan=(0.0,6.0)
# p=[σ=>10.0,ρ=>28.0,β=>8/3]
# prob=ODEProblem(sys,u,tspan,p)

# # Condition(u,t,integrator)=t-3
# # affect!(integrator)=integrator.u[4]=10
# # cb=ContinuousCallback(Condition,affect!)

# Condition1(u,t,integrator)=t-3
# affect1!(integrator)=integrator.u[4]=10
# cb1=ContinuousCallback(Condition1,affect1!)
# Condition2(u,t,integrator)=u[1]-5
# affect2!(integrator)=integrator.u[2]=10
# cb2=ContinuousCallback(Condition2,affect2!)
# cb=CallbackSet(cb1,cb2)

# sol=solve(prob,Tsit5(),callback=cb)
# println(sol)
# plot(sol)



#—————————————————————————————————————————— 离散回调 ————————————————————————————————————————
# using ModelingToolkit
# using DifferentialEquations
# using Plots

# @parameters σ,ρ,β
# @variables t,x(t),y(t),z(t)

# D=Differential(t)
# eqs=[D(x)~σ*(y-x),
#     D(y)~x*(ρ-z)-y,
#     D(z)~x*y-β*z]
# @named sys=ODESystem(eqs,t,[x,y,z],[σ,ρ,β])

# u=[x=>1.0,y=>0.0,z=>0]
# tspan=(0.0,6.0)
# p=[σ=>10,ρ=>28,β=>8/3]
# prob=ODEProblem(sys,u,tspan,p)

# # Condition(u,t,integrator)= t==3
# # affect!(integrator)=integrator.u[3]=10
# # cb=DiscreteCallback(Condition,affect!)
# # sol=solve(prob,Tsit5(),callback=cb,tstops=[3])
# time=[3,5]
# affect!(integrator)=integrator.u[2]=10
# cb=PresetTimeCallback(time,affect!)
# sol=solve(prob,Tsit5(),callback=cb)

# println(sol)
# plot(sol)



#—————————————————————————————————————————— 非线性方程 ——————————————————————————————————————————
# using ModelingToolkit
# using NonlinearSolve

# @variables x,y,z
# @parameters σ,ρ,β
# eqs=[0~σ*(y-x),
#      0~x*(ρ-z)-y,
#      0~x*y-β*z]
# ns=NonlinearSystem(eqs,[x,y,z],[σ,ρ,β])

# guess=[x=>1.0,y=>0.0,z=>0.0]
# ps=[σ=>10.0,ρ=>26.0,β=>8/3]
# prob=NonlinearProblem(ns,guess,ps)
# sol=solve(prob,NewtonRaphson())
# println(sol)



#—————————————————————————————————————————— 物性调用 ——————————————————————————————————————————
# using PyCall
# @pyimport CoolProp.CoolProp as CP
# rho=CP.PropsSI("D","P",101000,"Q",0,"Water")
# Temp=CP.PropsSI("T","P",101000,"Q",0,"Water")
# println([rho,Temp])

using ModelingToolkit, DifferentialEquations
using PyCall
@pyimport CoolProp.CoolProp as CP
@variables t T(t)
D = Differential(t)
function Prop(a,b,c,output) 
    if c==0
        medium="steel"
    elseif c==1
        medium="Water"
    elseif c==2
        medium="Air"
    elseif c==3
        medium="O2"
    elseif c==4
        medium="N2" 
    elseif c==5
        medium="H2"
    else
        println("error")
    end
    if c==0
        cp=350
    else
        if a==0 || a==1
            rho=CP.PropsSI("D","T",b,"Q",a,medium)
            cv=CP.PropsSI("O","T",b,"Q",a,medium)
            cp=CP.PropsSI("C","T",b,"Q",a,medium)
            thermalconductivity=CP.PropsSI("L","T",b,"Q",a,medium)
            visosity=CP.PropsSI("V","T",b,"Q",a,medium) 
            H=CP.PropsSI("Helmholtzmass","T",b,"Q",a,medium)
            p=CP.PropsSI("P","T",a,"Q",b,medium)
        else
            rho=CP.PropsSI("D","P",a,"T",b,medium)
            cv=CP.PropsSI("O","P",a,"T",b,medium)
            cp=CP.PropsSI("C","P",a,"T",b,medium)
            thermalconductivity=CP.PropsSI("L","P",a,"T",b,medium)
            visosity=CP.PropsSI("V","P",a,"T",b,medium) 
            H=CP.PropsSI("Helmholtzmass","P",a,"T",b,medium)
        end
    end
    if c==0
        M=0.056
    elseif c==1
        M=0.018
    elseif c==2
        M=0.032*0.79+0.028*0.21 
    elseif c==3
        M=0.032
    elseif c==4
        M=0.028 
    elseif c==5
        M=0.002
    else
        println("error")
    end   
    Rg=8.3145/M

    if output==1
        return(rho)
    elseif output==2
        return(cv)
    elseif output==3
        return(cp)
    elseif output==4
        return(thermalconductivity)
    elseif output==5
        return(visosity)
    elseif output==6
        return(H)
    elseif output==7
        return(p)
    elseif output==8
        return(Rg)
    else
        println("error")
    end          
end
@register Prop(a,b,c,switch)
eqs = [D(T) ~ Prop(t+300,1.0E5,4,8)]
@named sys = ODESystem(eqs)
sys = structural_simplify(sys)
u0=[T => 0]
tspan = (0.0,100.0)
prob = ODEProblem(sys,u0,tspan,[])
sol = solve(prob)
print(sol)
using Plots
plot(sol)



#———————————————————————————————————————— 饱和水蒸汽压 ————————————————————————————————————————
# function psat(T) "饱和水蒸汽压"
#     T=T-273.15
#     # a=-2.1794+0.02953*T-9.1837*10^-5*T^2+1.4454*10^-7*T^3
#     a=-2.18+0.0295*T-9.18*10^-5*T^2+1.44*10^-7*T^3
#     p=10^(a+5)
# end
# p=psat(298.15)
# println(p)



# ————————————————————————————————————————常微分方程优化 ————————————————————————————————————————
# using DifferentialEquations, DiffEqFlux, Plots
# using Plots

# function lotka_volterra!(du, u, p, t)
#   x, y = u
#   α, β, δ, γ = p
#   du[1] = dx = α*x - β*x*y
#   du[2] = dy = -δ*y + γ*x*y
# end

# u0 = [1.0, 1.0]
# tspan = (0.0, 10.0)
# tsteps = 0.0:0.1:10.0
# p = [1.5, 1.0, 3.0, 1.0]
# prob = ODEProblem(lotka_volterra!, u0, tspan, p)
# sol = solve(prob, Tsit5())
# plot(sol)

# function loss(p)
#   sol = solve(prob, Tsit5(), p=p, saveat = tsteps)
#   loss = sum(abs2, sol.-1)
#   return loss, sol
# end

# callback = function (p, l, pred)
#   display(l)
#   plt = plot(pred, ylim = (0, 6))
#   display(plt)
#   # Tell sciml_train to not halt the optimization. If return true, then
#   # optimization stops.
#   return false
# end

# result_ode = DiffEqFlux.sciml_train(loss, p, cb = callback)
# remade_solution = solve(remake(prob, p = result_ode.u), Tsit5(), saveat = tsteps)
# plot(remade_solution, ylim = (0, 6))



