using ReachabilityAnalysis
using ReachabilityBase.Arrays: SingleEntryVector

# parameters of the model
const g = 9.81            # gravity constant in m/s^2
const R = 0.1             # radius of center mass in m
const l = 0.5             # distance of rotors to center mass in m
const Mrotor = 0.1        # rotor mass in kg
const M = 1.0             # center mass in kg
const m = M + 4 * Mrotor  # total mass in kg
const mg = m * g;

# moments of inertia
const Jx = 0.4 * M * R^2 + 2 * l^2 * Mrotor
const Jy = Jx
const Jz = 0.4 * M * R^2 + 4 * l^2 * Mrotor
const Cyzx = (Jy - Jz) / Jx
const Czxy = (Jz - Jx) / Jy
const Cxyz = 0.0; #(Jx - Jy)/Jz

# control parameters
const u₁ = 1.0
const u₂ = 0.0
const u₃ = 0.0

@taylorize function quadrotor!(dx, x, p, t)
    x₁, x₂, x₃, x₄, x₅, x₆, x₇, x₈, x₉, x₁₀, x₁₁, x₁₂ = x

    # equations of the controllers
    F = (mg - 10 * (x₃ - u₁)) + 3 * x₆  # height control
    τϕ = -(x₇ - u₂) - x₁₀               # roll control
    τθ = -(x₈ - u₃) - x₁₁               # pitch control
    local τψ = 0.0                    # heading is uncontrolled

    Tx = τϕ / Jx
    Ty = τθ / Jy
    Tz = τψ / Jz
    F_m = F / m

    # some abbreviations
    sx7 = sin(x₇)
    cx7 = cos(x₇)
    sx8 = sin(x₈)
    cx8 = cos(x₈)
    sx9 = sin(x₉)
    cx9 = cos(x₉)

    sx7sx9 = sx7 * sx9
    sx7cx9 = sx7 * cx9
    cx7sx9 = cx7 * sx9
    cx7cx9 = cx7 * cx9
    sx7cx8 = sx7 * cx8
    cx7cx8 = cx7 * cx8
    sx7_cx8 = sx7 / cx8
    cx7_cx8 = cx7 / cx8

    x4cx8 = cx8 * x₄

    p11 = sx7_cx8 * x₁₁
    p12 = cx7_cx8 * x₁₂
    xdot9 = p11 + p12

    # differential equations for the quadrotor
    dx[1] = (cx9 * x4cx8 + (sx7cx9 * sx8 - cx7sx9) * x₅) + (cx7cx9 * sx8 + sx7sx9) * x₆
    dx[2] = (sx9 * x4cx8 + (sx7sx9 * sx8 + cx7cx9) * x₅) + (cx7sx9 * sx8 - sx7cx9) * x₆
    dx[3] = (sx8 * x₄ - sx7cx8 * x₅) - cx7cx8 * x₆
    dx[4] = (x₁₂ * x₅ - x₁₁ * x₆) - g * sx8
    dx[5] = (x₁₀ * x₆ - x₁₂ * x₄) + g * sx7cx8
    dx[6] = (x₁₁ * x₄ - x₁₀ * x₅) + (g * cx7cx8 - F_m)
    dx[7] = x₁₀ + sx8 * xdot9
    dx[8] = cx7 * x₁₁ - sx7 * x₁₂
    dx[9] = xdot9
    dx[10] = Cyzx * (x₁₁ * x₁₂) + Tx
    dx[11] = Czxy * (x₁₀ * x₁₂) + Ty
    dx[12] = Cxyz * (x₁₀ * x₁₁) + Tz
    return dx
end;

const T = 5.0
const v3 = SingleEntryVector(3, 12, 1.0)

@inline function quad_property(sol)
    tf = tend(sol)

    # Condition: b1 = (x[3] < 1.4) for all time
    unsafe1 = HalfSpace(-v3, -1.4) # unsafe: -x3 <= -1.4
    b1 = all([isdisjoint(unsafe1, set(R)) for R in sol((0.0, tf))])
    #b1 = ρ(v3, sol) < 1.4

    # Condition: x[3] > 0.9 for t ≥ 1.0
    unsafe2 = HalfSpace(v3, 0.9) # unsafe: x3 <= 0.9
    b2 = all([isdisjoint(unsafe2, set(R)) for R in sol((1.0, tf))])

    # Condition: x[3] ⊆ Interval(0.98, 1.02) for t = 5.0
    b3 = set(project(sol[end]; vars=(3))) ⊆ Interval(0.98, 1.02)

    return b1 && b2 && b3
end

function quadrotor(; Wpos, Wvel)
    # initial condition
    ΔX0 = [Wpos, Wpos, Wpos, Wvel, Wvel, Wvel, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    X0 = Hyperrectangle(zeros(12), ΔX0)

    # initial-value problem
    prob = @ivp(x' = quadrotor!(x), dim:12, x(0) ∈ X0)

    return prob
end;

Wpos = 0.1
Wvel = 0.1
prob = quadrotor(; Wpos=Wpos, Wvel=Wvel)
alg = TMJets(; abstol=1e-7, orderT=5, orderQ=1, adaptive=false)

sol = solve(prob; T=T, alg=alg)
solz1 = overapproximate(sol, Zonotope);

@assert quad_property(solz1) "the property should be proven"

Wpos = 0.4
Wvel = 0.4
prob = quadrotor(; Wpos=Wpos, Wvel=Wvel)
alg = TMJets(; abstol=1e-7, orderT=5, orderQ=1, adaptive=false)

sol = solve(prob; T=T, alg=alg)
solz2 = overapproximate(sol, Zonotope);

@assert quad_property(solz2) "the property should be proven"

Wpos = 0.8
Wvel = 0.8
prob = quadrotor(; Wpos=Wpos, Wvel=Wvel)
alg = TMJets(; abstol=1e-7, orderT=5, orderQ=1, adaptive=false)

sol = solve(prob; T=T, alg=alg)
solz3 = overapproximate(sol, Zonotope);

@assert !quad_property(solz3) "the property should not be proven"
