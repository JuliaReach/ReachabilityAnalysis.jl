using ReachabilityAnalysis

const μ = 3.986e14 * 60^2
const r = 42164.0e3
const r² = r^2
const mc = 500.0
const n² = μ / r^3
const n = sqrt(n²)

const two_n = 2*n
const μ_r² = μ/r²;

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

    # x' = vx
    du[1] = vx

    # y' = vy
    du[2] = vy

    # vx' = n²x + 2n*vy + μ/(r^2) - μ/(rc^3)*(r+x) + ux/mc
    du[3] = (n²*x + two_n*vy) + ((μ_r² - μ_rc³*rx) + uxy[1])

    # vy' = n²y - 2n*vx - μ/(rc^3)y + uy/mc
    du[4] = (n²*y - two_n*vx) - (μ_rc³*y - uxy[2])

    # t' = 1
    du[5] = one(x)

    return du
end

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

    # x' = vx
    du[1] = vx

    # y' = vy
    du[2] = vy

    # vx' = n²x + 2n*vy + μ/(r^2) - μ/(rc^3)*(r+x) + ux/mc
    du[3] = (n²*x + two_n*vy) + ((μ_r² - μ_rc³*rx) + uxy[1])

    # vy' = n²y - 2n*vx - μ/(rc^3)y + uy/mc
    du[4] = (n²*y - two_n*vx) - (μ_rc³*y - uxy[2])

    # t' = 1
    du[5] = one(x)

    return du
end

@taylorize function spacecraft_aborting!(du, u, p, t)
    x, y, vx, vy, t = u

    rx = r + x
    rx² = rx^2
    y² = y^2
    rc = sqrt(rx² + y²)
    rc³ = rc^3
    μ_rc³ = μ / rc³

    # x' = vx
    du[1] = vx

    # y' = vy
    du[2] = vy

    # vx' = n²x + 2n*vy + μ/(r^2) - μ/(rc^3)*(r+x)
    du[3] = (n²*x + two_n*vy) + (μ_r² - μ_rc³*rx)

    # vy' = n²y - 2n*vx - μ/(rc^3)y
    du[4] = (n²*y - two_n*vx) - μ_rc³*y

    # t' = 1
    du[5] = one(x)

    return du
end

using Symbolics

const var = @variables x y vx vy t

function spacecraft(; abort_time=(120.0, 150.0))

    n = 4 + 1  # number of variables
    t_abort_lower, t_abort_upper = abort_time

    automaton = LightAutomaton(3)

    # mode 1 "approaching"
    invariant = HalfSpace(x <= -100, var)
    approaching = @system(x' = spacecraft_approaching!(x), dim:5, x ∈ invariant)

    # mode 2 "rendezvous attempt"
    invariant = HalfSpace(x >= -100, var)
    attempt = @system(x' = spacecraft_attempt!(x), dim:5, x ∈ invariant)

    # mode 3 "aborting"
    invariant = Universe(n)
    aborting = @system(x' = spacecraft_aborting!(x), dim: 5, x ∈ invariant)

    # transition "approach" -> "attempt"
    add_transition!(automaton, 1, 2, 1)
    guard = HalfSpace(x >= -100, var)
    t1 = @map(x -> x, dim: n, x ∈ guard)

    # transition "approach" -> "abort"
    add_transition!(automaton, 1, 3, 2)
    guard_time = HPolyhedron([t >= t_abort_lower, t <= t_abort_upper], var)
    t2 = @map(x -> x, dim: n, x ∈ guard_time)

    # transition "attempt" -> "abort"
    add_transition!(automaton, 2, 3, 3)
    t3 = @map(x -> x, dim: n, x ∈ guard_time)

    H = HybridSystem(automaton=automaton,
                     modes=[approaching, attempt, aborting],
                     resetmaps=[t1, t2, t3])

     # initial condition in mode 1
     X0 = Hyperrectangle([-900., -400., 0., 0., 0.], [25., 25., 0., 0., 0.])
     init = [(1, X0)]

     return InitialValueProblem(H, init)
end

const tan30 = tand(30)

LineOfSightCone = HPolyhedron([x >= -100, y >= x*tan30, -y >= x*tan30], var)

Target = Hyperrectangle(zeros(2), [0.2, 0.2]);

function line_of_sight(sol)
    all_idx = findall(x -> x == 2, location.(sol))  # attempt
    for idx in all_idx
        verif = all(set(R) ⊆ LineOfSightCone for R in sol[idx][2:end])
        !verif && return false
    end
    return true
end

function absolute_velocity(R)
    vx, vy = 3, 4
    vx2 = set(overapproximate(project(R, vars=(vx,)), Interval))
    vy2 = set(overapproximate(project(R, vars=(vy,)), Interval))
    return sqrt(vx2.dat + vy2.dat)
end

function velocity_constraint(sol)
    all_idx = findall(x -> x == 2, location.(sol)) # attempt mode
    for idx in all_idx
        # maximum velocity measured in m/min
        verif = all(absolute_velocity(R) < 0.055 * 60. for R in sol[idx])
        !verif && return false
    end
    return true
end

function target_avoidance(sol)
    all_idx = findall(x -> x == 3, location.(sol)) # aborting mode
    for idx in all_idx
        verif = all(is_intersection_empty(set(Projection(R, [x, y])), Target) for R in sol[idx])
        !verif && return false
    end
    return true
end

using Plots

function solve_spacecraft(prob; k=25, s=missing)

    # template directions
    boxdirs = BoxDirections(5)

    # transition from mode 1 to mode 2
    sol12 = solve(prob,
                tspan=(0.0, 200.0),
                alg=TMJets(abstol=1e-5, maxsteps=10_000, orderT=5, orderQ=1, disjointness=BoxEnclosure()),
                max_jumps=1,
                intersect_source_invariant=false,
                intersection_method=TemplateHullIntersection(boxdirs),
                clustering_method=LazyClustering(1),
                disjointness_method=BoxEnclosure())

    sol12jump = overapproximate(sol12[2](120 .. 150), Zonotope)
    t0 = tstart(sol12jump[1])
    sol12jump_c = cluster(sol12jump, 1:length(sol12jump), BoxClustering(k, s))

    # transition from mode 2 to mode 3
    H = system(prob)
    sol3 = solve(IVP(mode(H, 3), [set(X) for X in sol12jump_c]),
                 tspan=(t0, 200.0),
                 alg=TMJets(abstol=1e-10, orderT=7, orderQ=1, disjointness=BoxEnclosure()))
    d = Dict{Symbol, Any}(:loc_id => 3)

    F12 = [fp for fp in sol12.F]
    F23 = [Flowpipe(fp.Xk, d) for fp in sol3.F]
    return HybridFlowpipe(vcat(F12, F23))
end

prob = spacecraft()
sol = solve_spacecraft(prob; k=25, s=missing)
solz = overapproximate(sol, Zonotope);

plot(solz, vars=(1, 2), xlab="x", ylab="y")

idx_approaching = findall(x -> x == 1, location.(solz))
idx_attempt = findall(x -> x == 2, location.(solz))
idx_aborting = findall(x -> x == 3, location.(solz))

using Plots, Plots.PlotMeasures, LaTeXStrings

fig = plot(legend=:bottomright, tickfont=font(10, "Times"), guidefontsize=15,
           xlab=L"x", ylab=L"y", lw=0.0,
           xtick=[-750, -500, -250, 0, 250.], ytick=[-400, -300, -200, -100, 0.],
           xlims=(-1000.0, 400.0), ylims=(-450.0, 0.0), size=(600, 600),
           bottom_margin=6mm, left_margin=2mm, right_margin=4mm, top_margin=3mm)

plot!(fig, solz[idx_approaching[1]], vars=(1, 2), lw=0.0, color=:lightgreen, alpha=1., lab="Approaching")
for k in 2:length(idx_approaching)
    plot!(fig, solz[idx_approaching[k]], vars=(1, 2), lw=0.0, color=:lightgreen, alpha=1.)
end

plot!(fig, solz[idx_attempt[1]], vars=(1, 2), lw=0.0, color=:red, alpha=1., lab="Attempt")
for k in 2:length(idx_attempt)
    plot!(fig, solz[idx_attempt[k]], vars=(1, 2), lw=0.0, color=:red, alpha=1.)
end

plot!(fig, solz[idx_aborting[1]], vars=(1, 2), lw=0.0, color=:cyan, alpha=1., lab="Aborting")
for k in 2:length(idx_aborting)
    plot!(fig, solz[idx_aborting[k]], vars=(1, 2), lw=0.0, color=:cyan, alpha=1.)
end

fig
