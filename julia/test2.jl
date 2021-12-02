using ModelingToolkit, DifferentialEquations
# Basic electric components
@parameters t

@connector function Pin(;name)
    @variables v(t) i(t)
    ODESystem(Equation[], t, [v, i], [], name=name, defaults=[v=>1.0, i=>1.0])
end

function ModelingToolkit.connect(::Type{Pin}, ps...)
    eqs = [
           0 ~ sum(p->p.i, ps) # KCL
          ]
    # KVL
    for i in 1:length(ps)-1
        push!(eqs, ps[i].v ~ ps[i+1].v)
    end

    return eqs
end

function Ground(;name)
    @named g = Pin()
    eqs = [g.v ~ 0]
    ODESystem(eqs, t, [], [], systems=[g], name=name)
end

function Resistor(;name, R = 1.0)
    val = R
    @named p = Pin()
    @named n = Pin()
    @variables v(t)
    @parameters R
    eqs = [
           v ~ p.v - n.v
           0 ~ p.i + n.i
           v ~ p.i * R 
          ]
    ODESystem(eqs, t, [v], [R], systems=[p, n], defaults=Dict(R => val), name=name)
end
function Capacitor(;name, C = 1.0)
    val = C
    @named p = Pin()
    @named n = Pin()
    @variables v(t)
    @parameters C
    D = Differential(t)
    eqs = [
           v ~ p.v - n.v
           0 ~ p.i + n.i
           D(v) ~ p.i / C
          ]
    ODESystem(eqs, t, [v], [C], systems=[p, n], defaults=Dict(C => val), name=name)
end

function Inductor(; name, L = 1.0)
    val = L
    @named p = Pin()
    @named n = Pin()
    @variables v(t) i(t)
    @parameters L
    D = Differential(t)
    eqs = [
           v ~ p.v - n.v
           0 ~ p.i + n.i
           i ~ p.i
           D(i) ~ v / L
          ]
    ODESystem(eqs, t, [v, i], [L], systems=[p, n], defaults=Dict(L => val), name=name)
end


function ChangeableVoltage(;name,vol=16.0)
    @named p = Pin()
    @named n = Pin()
    @variables v(t)
    @parameters Vol
    eqs = [
           v ~ Vol*sin(2π*t)
           v ~ p.v - n.v
           0 ~ p.i + n.i
          ]
    ODESystem(eqs, t, [v], [Vol], systems=[p, n], defaults=Dict(Vol => vol),name=name)
end

@named resistor = Resistor(R=3.0)
@named capacitor = Capacitor(C=1.0/24)
@named source = ChangeableVoltage(vol=10.0)
@named inductor = Inductor(L=0.1)
@named ground = Ground()


rc_eqs = [
          connect(source.p, capacitor.p)
          connect(capacitor.n, inductor.p)
          connect(inductor.n, resistor.p)
          connect(source.n,resistor.n,ground.g)
         ]

sys = [resistor, inductor, capacitor,source,ground]
@named rc_model = ODESystem(rc_eqs, t, systems=sys)
sys = structural_simplify(rc_model)

u0 = [
    capacitor.v => 0.0
    capacitor.p.i => 0.0
    inductor.i => 0
    inductor.v => 0
     ]
P=[3.0,1.0/24,0.1,10.0] #系统的参数
prob = ODAEProblem(sys, u0, (0, 10.0),P)
sol = solve(prob, Tsit5())

using Plots
p1 = plot(sol,vars=[capacitor.v capacitor.p.i],xlims = (0,10),ylim = (-10,15))
p2 = plot(sol,vars=[inductor.v inductor.i],xlims = (0,10),ylim = (-5,5))
plot(p1,p2,layout=(2,1))

parameters(sys)
states(sys)

using RecursiveArrayTools
t = collect(range(0,stop=10,length=1000))
randomized = VectorOfArray([(sol(t[i]) + .5randn(2)) for i in 1:length(t)])
data = convert(Array,randomized) 
p1=plot(t,data[1,:],ylim=(-8,8))
p1=plot!(t,data[2,:],ylim=(-8,8))
p2=plot(sol,vars=[inductor.i,capacitor.v],ylim=(-8,8))
plot(p1,p2,layout=(2,1))

using DiffEqParamEstim
export DiffEqObjective, DiffEqParamEstim, L2Loss, LogLikeLoss, Regularization, TwoStageCost, build_loss_objective, build_lsoptim_objective, colloc_grad, l2lossgradient!, lm_fit, multiple_shooting_objective, prior_loss, two_stage_method
cost_function = build_loss_objective(prob,Tsit5(),L2Loss(t,data),
                                     maxiters=10000,verbose=false)
vals = 0:0.1:20.0
plot(vals,[cost_function([3.0,1.0/24,0.1,i]) for i in vals],yscale=:log10,
     xaxis = "Parameter", yaxis = "Cost", title = "Cost Function",
     lw = 3)

using Optim
result = optimize(cost_function, [2.0,0.1,0.5,16])
result.minimizer
