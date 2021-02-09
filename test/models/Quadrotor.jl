using ReachabilityAnalysis
using ReachabilityAnalysis: is_intersection_empty

#The task is to change the height from ``0``~[m] to ``1``~[m] within ``5``~[s].

const g = 9.81           # gravity constant in m/s^2
const R = 0.1            # radius of center mass in m
const l = 0.5            # distance of motors to center mass in m
const Mrotor = 0.1       # motor mass in kg
const M = 1.0            # center mass in kg
const m = M + 4*Mrotor   # total mass in kg
const mg = m*g

const Jx = (2/5)*M*R^2 + 2*l^2*Mrotor
const Jy = Jx
const Jz = (2/5)*M*R^2 + 4*l^2*Mrotor
const Cyzx = (Jy - Jz)/Jx
const Czxy = (Jz - Jx)/Jy
const Cxyz = 0.0 #(Jx - Jy)/Jz

const u₁ = 1.0
const u₂ = 0.0
const u₃ = 0.0

const Tspan = (0.0, 5.0)
const v3 = LazySets.SingleEntryVector(3, 12, 1.0)

@inline function quad_property(solz)
    tf = tend(solz)

    #Condition: b1 = (x[3] < 1.4) for all time
    unsafe1 = HalfSpace(-v3, -1.4) # unsafe: -x3 <= -1.4
    b1 =  all([is_intersection_empty(unsafe1, set(R)) for R in solz(0.0 .. tf)])
    #b1 = ρ(v3, solz) < 1.4

    #Condition: x[3] > 0.9 for t ≥ 1.0
    unsafe2 = HalfSpace(v3, 0.9) # unsafe: x3 <= 0.9
    b2 = all([is_intersection_empty(unsafe2, set(R)) for R in solz(1.0 .. tf)])

    #Condition: x[3] ⊆ Interval(0.98, 1.02) for t ≥ 5.0 (t=5 is `tf`)
    b3 = set(project(solz[end], vars=(3))) ⊆ Interval(0.98, 1.02)

    return b1 && b2 && b3
end

@taylorize function quadrotor!(dx, x, params, t)
    #unwrap the variables and the controllers; the last three are the controllers
    #x₁, x₂, x₃, x₄, x₅, x₆, x₇, x₈, x₉, x₁₀, x₁₁, x₁₂, u₁, u₂, u₃ = x
    x₁  = x[1]
    x₂  = x[2]
    x₃  = x[3]
    x₄  = x[4]
    x₅  = x[5]
    x₆  = x[6]
    x₇  = x[7]
    x₈  = x[8]
    x₉  = x[9]
    x₁₀ = x[10]
    x₁₁ = x[11]
    x₁₂ = x[12]

    #equations of the controllers
    F = (mg - 10*(x₃ - u₁)) + 3*x₆  # height control
    τϕ = -(x₇ - u₂) - x₁₀            # roll control
    τθ = -(x₈ - u₃) - x₁₁            # pitch control
    local τψ = 0.0                   # heading is uncontrolled

    Tx = τϕ/Jx
    Ty = τθ/Jy
    Tz = τψ/Jz
    F_m = F/m

    #Some abbreviations
    sx7 = sin(x₇)
    cx7 = cos(x₇)
    sx8 = sin(x₈)
    cx8 = cos(x₈)
    sx9 = sin(x₉)
    cx9 = cos(x₉)

    sx7sx9 = sx7*sx9
    sx7cx9 = sx7*cx9
    cx7sx9 = cx7*sx9
    cx7cx9 = cx7*cx9
    sx7cx8 = sx7*cx8
    cx7cx8 = cx7*cx8
    sx7_cx8 = sx7/cx8
    cx7_cx8 = cx7/cx8

    x4cx8 = cx8*x₄

    p11 = sx7_cx8*x₁₁
    p12 = cx7_cx8*x₁₂
    xdot9 = p11 + p12

    #differential equations for the quadrotor

    dx[1] = (cx9*x4cx8 + (sx7cx9*sx8 - cx7sx9)*x₅) + (cx7cx9*sx8 + sx7sx9)*x₆
    dx[2] = (sx9*x4cx8 + (sx7sx9*sx8 + cx7cx9)*x₅) + (cx7sx9*sx8 - sx7cx9)*x₆
    dx[3] = (sx8*x₄ - sx7cx8*x₅) - cx7cx8*x₆
    dx[4] = (x₁₂*x₅ - x₁₁*x₆) - g*sx8
    dx[5] = (x₁₀*x₆ - x₁₂*x₄) + g*sx7cx8
    dx[6] = (x₁₁*x₄ - x₁₀*x₅) + (g*cx7cx8 - F_m)
    dx[7] = x₁₀ + sx8*xdot9
    dx[8] = cx7*x₁₁ - sx7*x₁₂
    dx[9] = xdot9
    dx[10] = Cyzx * (x₁₁ * x₁₂) + Tx
    dx[11] = Czxy * (x₁₀ * x₁₂) + Ty
    dx[12] = Cxyz * (x₁₀ * x₁₁) + Tz

    return dx
end

function quadrotor(; T=5.0, plot_vars=[0, 3],
                property=quad_property,
                project_reachset=true,
                Wpos = 0.4, Wvel = 0.4)

    #initial conditions
    X0c = zeros(12)
    ΔX0 = [Wpos, Wpos, Wpos, Wvel, Wvel, Wvel, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    X0 = Hyperrectangle(X0c, ΔX0)

    #initial-value problem
    prob = @ivp(x' = quadrotor!(x), dim: 12, x(0) ∈ X0);

    return prob
end

cases = ["Δ=0.1", "Δ=0.4", "Δ=0.8"];

Wpos = 0.1
Wvel = 0.1
prob = quadrotor(project_reachset=false, Wpos=Wpos, Wvel=Wvel)
alg = TMJets(abs_tol=1e-7, orderT=5, orderQ=1, adaptive=false)

sol1 = solve(prob, tspan=Tspan, alg=alg);
solz1 = overapproximate(sol1, Zonotope);

Wpos = 0.4
Wvel = 0.4
prob = quadrotor(project_reachset=false, Wpos=Wpos, Wvel=Wvel)
alg = TMJets(abs_tol=1e-7, orderT=5, orderQ=1, adaptive=false)

sol2 = solve(prob, tspan=Tspan, alg=alg);
solz2 = overapproximate(sol2, Zonotope);

#property = quad_property(solz2)
#println("Validate property, case $(cases[2]) : $(property)")

Wpos = 0.8
Wvel = 0.8
prob = quadrotor(project_reachset=false, Wpos=Wpos, Wvel=Wvel)
alg = TMJets(abs_tol=1e-7, orderT=5, orderQ=1, adaptive=false)

sol3 = solve(prob, tspan=Tspan, alg=alg);
solz3 = overapproximate(sol3, Zonotope);

using Plots

Plots.plot(solz3,  vars=(0, 3), linecolor="green",  color=:green,  alpha=0.8)
Plots.plot!(solz2, vars=(0, 3), linecolor="blue",   color=:blue,   alpha=0.8)
Plots.plot!(solz1, vars=(0, 3), linecolor="yellow", color=:yellow, alpha=0.8,
    xlab="t", ylab="x3",
    xtick=[0., 1., 2., 3., 4., 5.], ytick=[-1., -0.5, 0., 0.5, 1., 1.5],
    xlims=(0., 5.), ylims=(-1., 1.5))

