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



# 常微分方程
# using ModelingToolkit
# using OrdinaryDiffEq
# using Plots
# @parameters σ,ρ,β
# @variables t,x(t),y(t),z(t)
# D=Differential(t)
# eqs=(D(x)~σ*(y-x),
#     D(y)~x*(ρ-z)-y,
#     D(z)~x*y-β*z)
# @named sys=ODESystem(eqs,t,[x,y,z],[σ,ρ,β])

# u=[1.0,0.0,0.0]
# tspan=(0.0,6.0)
# p=[10.0,28.0,8/3]
# prob=ODEProblem(sys,u,tspan,p)
# sol=solve(prob,Tsit5())
# println(sol)
# plot(sol)



#微分代数方程-中间变量
# using ModelingToolkit
# using OrdinaryDiffEq
# using Plots

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



# 微分代数方程-外部函数
# using ModelingToolkit
# using OrdinaryDiffEq
# using Plots

# @parameters σ,ρ,β
# @variables t,x(t),y(t),z(t),f(t)
# D=Differential(t)
# eqs=[f~10+t/1000,
#     D(x)~f+σ*(y-x),
#     D(y)~x*(ρ-z)-y,
#     D(z)~x*y-β*z]
# @named sys=ODESystem(eqs)
# sys_simple=structural_simplify(sys)

# u=[x=>1.0,y=>0.0,z=>0.0]
# tspan=(0.0,6.0)
# p=[σ=>10.0,ρ=>28.0,β=>8/3]
# prob=ODEProblem(sys_simple,u,tspan,p)
# sol=solve(prob,Tsit5())
# println(sol)
# println(sol[f])
# plot(sol,vars=[x,y,z,f])



# 连续回调
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

using ModelingToolkit
using DifferentialEquations
using PyCall
using Plots

@parameters σ,ρ,β
@variables t,x(t),y(t),z(t)
@pyimport CoolProp.CoolProp as CP
D=Differential(t)
eqs=(D(x)~ω/1000,
    D(y)~1000,
    D(z)~0)
@named sys=ODESystem(eqs,t,[x,y,z],[ω])

u=[x=>1.0,y=>101000,z=>0.0]
tspan=(0.0,6.0)
p=[ω=>CP.PropsSI("D","P",101000,"Q",0,"Water")]
prob=ODEProblem(sys,u,tspan,p)

Condition(u,t,integrator)=t
affect!(integrator)=integrator.p[1]=CP.PropsSI("D","P",y,"Q",0,"Water")+1
cb=ContinuousCallback(Condition,affect!)

# Condition1(u,t,integrator)=t-3
# affect1!(integrator)=integrator.u[4]=10
# cb1=ContinuousCallback(Condition1,affect1!)
# Condition2(u,t,integrator)=u[1]-5
# affect2!(integrator)=integrator.u[2]=10
# cb2=ContinuousCallback(Condition2,affect2!)
# cb=CallbackSet(cb1,cb2)

sol=solve(prob,Tsit5(),callback=cb)
println(sol[ω])
plot(sol)



# 离散回调
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



# 非线性方程
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



# 物性调用
# using PyCall
# @pyimport CoolProp.CoolProp as CP
# rho=CP.PropsSI("D","P",101000,"Q",0,"Water")
# Temp=CP.PropsSI("T","P",101000,"Q",0,"Water")
# println([rho,Temp])



# function psat(T) "饱和水蒸汽压"
#     T=T-273.15
#     # a=-2.1794+0.02953*T-9.1837*10^-5*T^2+1.4454*10^-7*T^3
#     a=-2.18+0.0295*T-9.18*10^-5*T^2+1.44*10^-7*T^3
#     p=10^(a+5)
# end
# p=psat(298.15)
# println(p)

