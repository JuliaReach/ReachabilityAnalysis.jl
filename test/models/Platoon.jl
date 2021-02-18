using ReachabilityAnalysis, SparseArrays, ModelingToolkit

const var = @variables x[1:9] t

function platoon_connected(; deterministic_switching::Bool=true,
                             c1=5.0)  # clock constraints
    n = 9 + 1

    # x' = Ax + Bu + c
    A = Matrix{Float64}(undef, n, n)
    A[1, :] = [0, 1.0, 0, 0, 0, 0, 0, 0, 0, 0]
    A[2, :] = [0, 0, -1.0, 0, 0, 0, 0, 0, 0, 0]
    A[3, :] = [1.6050, 4.8680, -3.5754, -0.8198, 0.4270, -0.0450, -0.1942,  0.3626, -0.0946, 0.]
    A[4, :] = [0, 0, 0, 0, 1.0, 0, 0, 0, 0, 0,]
    A[5, :] = [0, 0, 1.0, 0, 0, -1.0, 0, 0, 0, 0]
    A[6, :] = [0.8718, 3.8140, -0.0754,  1.1936, 3.6258, -3.2396, -0.5950,  0.1294, -0.0796, 0.]
    A[7, :] = [0, 0, 0, 0, 0, 0, 0, 1.0, 0, 0]
    A[8, :] = [0, 0, 0, 0, 0, 1.0, 0, 0, -1.0, 0]
    A[9, :] = [0.7132, 3.5730, -0.0964,  0.8472, 3.2568, -0.0876,  1.2726,  3.0720, -3.1356, 0.]
    A[10, :] = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0]; # t' = 1

    if deterministic_switching
        invariant = HalfSpace(t <= c1, var)
    else
        invariant = Universe(n)
    end

    # acceleration of the lead vehicle + time
    B = sparse([2], [1], [1.0], n, 1)
    U = Hyperrectangle(low=[-9.], high=[1.])
    c = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1.0]
    @system(x' = Ax + Bu + c, x ∈ invariant, u ∈ U)
end

function platoon_disconnected(; deterministic_switching::Bool=true,
                                c2=5.0)  # clock constraints
    n = 10 # 9 dimensions + time

    # x' = Ax + Bu + c
    A = Matrix{Float64}(undef, n, n)
    A[1, :] = [0, 1.0, 0, 0, 0, 0, 0, 0, 0, 0]
    A[2, :] = [0, 0, -1.0, 0, 0, 0, 0, 0, 0, 0]
    A[3, :] = [1.6050, 4.8680, -3.5754, 0, 0, 0, 0, 0, 0, 0]
    A[4, :] = [0, 0, 0, 0, 1.0, 0, 0, 0, 0, 0,]
    A[5, :] = [0, 0, 1.0, 0, 0, -1.0, 0, 0, 0, 0]
    A[6, :] = [0, 0, 0,  1.1936, 3.6258, -3.2396, 0, 0, 0, 0.]
    A[7, :] = [0, 0, 0, 0, 0, 0, 0, 1.0, 0, 0]
    A[8, :] = [0, 0, 0, 0, 0, 1.0, 0, 0, -1.0, 0]
    A[9, :] = [0.7132, 3.5730, -0.0964,  0.8472, 3.2568, -0.0876,  1.2726,  3.0720, -3.1356, 0.]
    A[10, :] = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0]; # t' = 1

    if deterministic_switching
        invariant = HalfSpace(t <= c2, var)
    else
        invariant = Universe(n)
    end

    # acceleration of the lead vehicle + time
    B = sparse([2], [1], [1.0], n, 1)
    U = Hyperrectangle(low=[-9.], high=[1.])
    c = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1.0]
    @system(x' = Ax + Bu + c, x ∈ invariant, u ∈ U)
end

function platoon(; deterministic_switching::Bool=true,
                   c1=5.0,  # clock constraints
                   c2=5.0,  # clock constraints
                   tb=10.0,  # lower bound for loss of communication
                   tc=20.0, tr=20.0) # upper bound for loss of communication (tc) and reset time (tr)

    # three variables for each vehicle, (ei, d(et)/dt, ai) for
    # (spacing error, relative velocity, speed), and the last dimension is time
    n = 9 + 1

    # transition graph
    automaton = LightAutomaton(2)
    add_transition!(automaton, 1, 2, 1)
    add_transition!(automaton, 2, 1, 2)

    # modes
    mode1 = platoon_connected(deterministic_switching=deterministic_switching, c1=c1)
    mode2 = platoon_disconnected(deterministic_switching=deterministic_switching, c2=c2)
    modes = [mode1, mode2]

    # common reset
    reset = Dict(n => 0.)

    # transition l1 -> l2
    if deterministic_switching
        guard = Hyperplane(t == c1, var)
    else
        guard = HPolyhedron([tb <= t, t <= tc], var)
    end
    t1 = ConstrainedResetMap(n, guard, reset)

    # transition l2 -> l1
    if deterministic_switching
        guard = Hyperplane(t == c2, var)
    else
        guard = HalfSpace(t <= tr, var)
    end
    t2 = ConstrainedResetMap(n, guard, reset)
    resetmaps = [t1, t2]

    H = HybridSystem(automaton, modes, resetmaps, [AutonomousSwitching()])

    # initial condition is at the orgin in mode 1
    X0 = BallInf(zeros(n), 0.0)
    initial_condition = [(1, X0)]

    return IVP(H, initial_condition)
end

function dmin_specification(sol, dmin)
    return (-ρ(sparsevec([1], [-1.0], 10), sol) > -dmin) &&
           (-ρ(sparsevec([4], [-1.0], 10), sol) > -dmin) &&
           (-ρ(sparsevec([7], [-1.0], 10), sol) > -dmin)
end

prob_PLAD01 = platoon();

const boxdirs = BoxDirections(10)
const octdirs = OctDirections(10);

length(boxdirs)

alg = BOX(δ=0.01)
sol_PLAD01_BND42 = solve(prob_PLAD01,
                         alg=alg,
                         clustering_method=BoxClustering(1),
                         intersection_method=TemplateHullIntersection(boxdirs),
                         intersect_source_invariant=false,
                         tspan = (0.0 .. 20.0));

dmin_specification(sol_PLAD01_BND42, 42)

-ρ(sparsevec([1], [-1.0], 10), sol_PLAD01_BND42)

-ρ(sparsevec([4], [-1.0], 10), sol_PLAD01_BND42)

-ρ(sparsevec([7], [-1.0], 10), sol_PLAD01_BND42)

using Plots, LaTeXStrings

plot(sol_PLAD01_BND42, vars=(0, 1), xlab=L"t", ylab=L"x_1", title="PLAD01 - BND42", lw=0.1)
plot!(x->x, x->-42., 0., 20., linewidth=2, color="red", linestyle=:dash, leg=nothing)

length(octdirs)


alg = LGG09(δ=0.03, template=octdirs, approx_model=Forward(setops=octdirs));
sol_PLAD01_BND30 = solve(prob_PLAD01,
                         alg=alg,
                         clustering_method=LazyClustering(1),
                         intersection_method=TemplateHullIntersection(octdirs),
                         intersect_source_invariant=false,
                         tspan = (0.0 .. 20.0));

dmin_specification(sol_PLAD01_BND30, 30)

-ρ(sparsevec([1], [-1.0], 10), sol_PLAD01_BND30)

-ρ(sparsevec([4], [-1.0], 10), sol_PLAD01_BND30)

-ρ(sparsevec([7], [-1.0], 10), sol_PLAD01_BND30)

plot(sol_PLAD01_BND30, vars=(0, 1), xlab=L"t", ylab=L"x_1", title="PLAD01 - BND30", lw=0.1)
plot!(x->x, x->-30., 0., 20., linewidth=2, color="red", linestyle=:dash, leg=nothing)

