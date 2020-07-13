# # Spacecraft
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/models/spacecraft.ipynb)
#
#md # !!! note "Overview"
#md #     System type: polynomial continuous system\
#md #     State dimension: 2\
#md #     Application domain: Chemical kinetics
#
# ## Model description
#
# A chemical reaction is said to be *autocatalytic* if one of the reaction products is
# also a catalyst for the same or a coupled reaction, and such a reaction is called an autocatalytic reaction.
# We refer to the wikipedia article [Autocatalysis](https://en.wikipedia.org/wiki/Autocatalysis) for details.

# The Brusselator is a mathematical model for a class of autocatalytic reactions.
# The dynamics of the Brusselator is given by the two-dimensional ODE
#
# ```math
#   \left\{ \begin{array}{lcl} \dot{x} & = & A + x^2\cdot y - B\cdot x - x \\
#    \dot{y} & = & B\cdot x - x^2\cdot y \end{array} \right.
# ```

# ## Reachability settings
#
# TODO:
# ISSC01 description
# The initial set ``X_0`` is defined by ``B = {x_i \in [-0001,+0.0001]``, ``\forall i}`` and
# ``U \in [0, 0.1] \times [0.8, 1.] \times [0.9, 1.]``. Then ``X_0 = B \times U``
# and ``u \in U``.
# These settings are taken from [1].
# ```math
#   \dot{x} = A \cdot x + B \cdot u
# ```


using ReachabilityAnalysis, SparseArrays
using ReachabilityAnalysis: TaylorModelReachSet, AbstractLazyReachSet

const μ = 3.986e14 * 60^2
const r = 42164.0e3
const r² = r^2
const mc = 500.0
const n² = μ / r^3
const n = sqrt(n²)

const two_n = 2*n
const μ_r² = μ/r²

# columns correspond to x, y, vx, vy, t
const K₁ = [-28.8287 0.1005 -1449.9754 0.0046 0.0;
            -0.087 -33.2562 0.00462 -1451.5013 0.0]
const K₂ = [-288.0288 0.1312 -9614.9898 0.0 0.0;
            -0.1312 -288.0 0.0 -9614.9883   0.0]

const K₁mc = K₁/mc
const K₂mc = K₂/mc

function mymul!(v, A, x)
    @inbounds for ind = 1:length(v)
        v[ind] = zero(x[1])
        for jind = 1:length(x)
            v[ind] += A[ind, jind] * x[jind]
        end
    end
    return nothing
end

# dynamics in the 'approaching' mode
@taylorize function spacecraft_approaching!(du, u, p, t)
    x, y, vx, vy, t = u

    rx = r + x
    rx² = rx^2
    y² = y^2
    rc = sqrt(rx² + y²)
    rc³ = rc^3
    μ_rc³ = μ / rc³

    uxy = Vector{typeof(x)}(undef, 2)
    mymul!(uxy, K₁mc, u)

    #x' = vx
    du[1] = vx

    #y' = vy
    du[2] = vy

    #vx' = n²x + 2n*vy + μ/(r^2) - μ/(rc^3)*(r+x) + ux/mc
    du[3] = (n²*x + two_n*vy) + ((μ_r² - μ_rc³*rx) + uxy[1])

    #vy' = n²y - 2n*vx - μ/(rc^3)y + uy/mc
    du[4] = (n²*y - two_n*vx) - (μ_rc³*y - uxy[2])

    #t' = 1
    du[5] = one(x)

    return du
end

# dynamics in the 'rendezvous attempt' mode
@taylorize function spacecraft_attempt!(du, u, p, t)
    x, y, vx, vy, t = u

    rx = r + x
    rx² = rx^2
    y² = y^2
    rc = sqrt(rx² + y²)
    rc³ = rc^3
    μ_rc³ = μ / rc³

    uxy = Vector{typeof(x)}(undef, 2)
    mymul!(uxy, K₂mc, u)

    #x' = vx
    du[1] = vx

    #y' = vy
    du[2] = vy

    #vx' = n²x + 2n*vy + μ/(r^2) - μ/(rc^3)*(r+x) + ux/mc
    du[3] = (n²*x + two_n*vy) + ((μ_r² - μ_rc³*rx) + uxy[1])

    #vy' = n²y - 2n*vx - μ/(rc^3)y + uy/mc
    du[4] = (n²*y - two_n*vx) - (μ_rc³*y - uxy[2])

    #t' = 1
    du[5] = one(x)

    return du
end

# dynamics in the 'aborting' mode
@taylorize function spacecraft_aborting!(du, u, p, t)
    x, y, vx, vy, t = u

    rx = r + x
    rx² = rx^2
    y² = y^2
    rc = sqrt(rx² + y²)
    rc³ = rc^3
    μ_rc³ = μ / rc³

    #x' = vx
    du[1] = vx

    #y' = vy
    du[2] = vy

    #vx' = n²x + 2n*vy + μ/(r^2) - μ/(rc^3)*(r+x)
    du[3] = (n²*x + two_n*vy) + (μ_r² - μ_rc³*rx)

    #vy' = n²y - 2n*vx - μ/(rc^3)y
    du[4] = (n²*y - two_n*vx) - μ_rc³*y

    #t' = 1
    du[5] = one(x)

    return du
end

function spacecraft(; X0 = Hyperrectangle([-900., -400., 0., 0., 0.],
                                          [25., 25., 0., 0., 0.]),
                         init=[(1, X0)],
                         abort_time=(120.0, 150.0))

    #variables
    x = 1   # x position
    y = 2   # y position
    vx = 3  # x velocity
    vy = 4  # y velocity
    t = 5   # time
    n = 4 + 1  # number of variables
    t_abort_lower, t_abort_upper = abort_time[1], abort_time[2]

    automaton = LightAutomaton(3)

    #mode 1 "approaching"
    invariant = HalfSpace(sparsevec([x], [1.], n), -100.) # x <= -100
    approaching = @system(x' = spacecraft_approaching!(x), dim:5, x ∈ invariant)

    #mode 2 ("rendezvous attempt")
    invariant = HalfSpace(sparsevec([x], [-1.], n), 100.) # x >= -100
    attempt = @system(x' = spacecraft_attempt!(x), dim:5, x ∈ invariant)

    #mode 3 "aborting"
    invariant = Universe(n)
    aborting = @system(x' = spacecraft_aborting!(x), dim: 5, x ∈ invariant)

    #transition "approach" -> "attempt"
    add_transition!(automaton, 1, 2, 1)
    guard = HalfSpace(sparsevec([x], [-1.], n), 100.) # x >= -100
    #t1 = @map(x -> x, dim: n, x ∈ guard)
    t1 = map(x -> x, guard)


    #transition "approach" -> "abort"
    add_transition!(automaton, 1, 3, 2)
    guard_time = HPolyhedron([HalfSpace(sparsevec([t], [-1.], n), -t_abort_lower),  # t >= t_abort_lower
                              HalfSpace(sparsevec([t], [1.], n), t_abort_upper)])   # t <= t_abort_upper
    #t2 = @map(x -> x, dim: n, x ∈ guard_time)
    t2 = map(x -> x, guard_time)

    #transition "attempt" -> "abort"
    add_transition!(automaton, 2, 3, 3)
    #t3 = @map(x -> x, dim: n, x ∈ guard_time)
    t3 = map(x -> x, guard_time)

    H = HybridSystem(automaton=automaton,
                     modes=[approaching, attempt, aborting],
                     resetmaps=[t1, t2, t3])

    return InitialValueProblem(H, init)
end



#-

# ## Results

# TODO: We use `LGG09` algorithm with sixth-order expansion in time and second order expansion
# in the spatial variables.

boxdirs = BoxDirections(5)

function solve_spacecraft(prob; k=25, s=missing)

    #transition from mode 1 to mode 2
    sol12 = solve(prob,
                tspan=(0.0, 200.0),
                alg=TMJets(abs_tol=1e-5, max_steps=10_000, orderT=5, orderQ=1, disjointness=BoxEnclosure()),
                max_jumps=1,
                intersect_source_invariant=false,
                intersection_method=TemplateHullIntersection(boxdirs),
                clustering_method=LazyClustering(1),
                disjointness_method=BoxEnclosure())

    sol12jump = overapproximate(sol12[2](120 .. 150), Zonotope)
    t0 = tstart(sol12jump[1])
    sol12jump_c = cluster(sol12jump, 1:length(sol12jump), BoxClustering(k, s))

    #transition from mode 2 to mode 3
    H = system(prob)
    sol3 = solve(IVP(mode(H, 3), [set(X) for X in sol12jump_c]),
                 tspan=(t0, 200.0),
                 alg=TMJets(abs_tol=1e-10, orderT=7, orderQ=1, disjointness=BoxEnclosure()))
    d = Dict{Symbol, Any}(:loc_id => 3)

    return HybridFlowpipe(vcat([fp for fp in sol12.F],
                               [Flowpipe(fp.Xk, d) for fp in sol3.F]))
end

prob = spacecraft()
sol = solve_spacecraft(prob; k=25, s=missing)
solz = overapproximate(sol, Zonotope)

#=
# verify that specifications hold
prop1 = line_of_sight(solz)
println("Line of sight property: $prop1")

prop2 = velocity_constraint(solz)
println("Velocity constraint property $prop2")

prop3 = target_avoidance(solz)
println("Target avoidance property: $prop3")

property = prop1 && prop2 && prop3
=#

#-

using Plots

idx_approaching = findall(x -> x == 1, location.(solz))
idx_attempt = findall(x -> x == 2, location.(solz))
idx_aborting = findall(x -> x == 3, location.(solz))

Plots.plot(legend=:bottomright)

for idx in idx_approaching
    Plots.plot!(solz[idx], vars=(1, 2), lw=0.0, color=:blue, alpha=1.,
    xlab="x", ylab="y",)
end
for idx in idx_attempt
    Plots.plot!(solz[idx], vars=(1, 2), lw=0.0, color=:red, alpha=1.)
end
for idx in idx_aborting
    Plots.plot!(solz[idx], vars=(1, 2), lw=0.0, color=:green, alpha=1.)
end
