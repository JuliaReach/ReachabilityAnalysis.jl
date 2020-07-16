# # Platooning
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/models/platoon.ipynb)
#
#md # !!! note "Overview"
#md #     System type: affine system\
#md #     State dimension: 9 + 1\
#md #     Application domain: Chemical kinetics
#
# ## Model description
#
# The platooning benchmark considers a platoon of three vehicles following each
# other. This benchmark considers loss of communication between vehicles. The
# initial discrete state is ``q_c``. Three scenarios are considered for the loss
# of communication:
#
#
# PLAA01 (arbitrary loss) The loss of communication can occur at any time. This
#       includes the possibility of no communication at all.
#
# PLADxy (loss at deterministic times) The loss of communication occurs at fixed
#       points in time, which are determined by clock constraints ``c_1`` and ``c_2``.
#       The clock t is reset when communication is lost and when it is re-established.
#       Note that the transitions have must-semantics, i.e., they take place as soon
#       as possible.
#
#       PLAD01: ``c_1 = c_2 = 5``.
#
# PLANxy (loss at nondeterministic times) The loss of communication occurs at
#       any time ``t ∈ [t_b, t_c]``. The clock t is reset when communication is lost
#       and when it is reestablished. Communication is reestablished at any time
#       ``t ∈ [0, t_r]``. This scenario covers loss of communication after an
#       arbitrarily long time ``t ≥ t_c`` by reestablishing communication in zero time.
#
#       PLAN01: ``t_b = 10, t_c = 20, t_r = 20``.

using ReachabilityAnalysis, SparseArrays, ModelingToolkit

# ## Dynamics of the "connected" platoon

function platoon_connected(; deterministic_switching::Bool=true,
                             c1=5.0)  # clock constraints
    n = 10 # 9 dimensions + time
    #x' = Ax + Bu + c
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
        invariant = HalfSpace(sparsevec([n], [1.], n), c1) # t <= c1
    else
        invariant = Universe(n)
    end

    #acceleration of the lead vehicle + time
    B = sparse([2], [1], [1.0], n, 1)
    U = Hyperrectangle(low=[-9.], high=[1.])
    c = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1.0]
    @system(x' = Ax + Bu + c, x ∈ invariant, u ∈ U)
end

# ## Dynamics of the "disconnected" platoon

function platoon_disconnected(; deterministic_switching::Bool=true,
                                c2=5.0)  # clock constraints
    n = 10 # 9 dimensions + time
    #x' = Ax + Bu + c
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
        invariant = HalfSpace(sparsevec([n], [1.], n), c2) # t <= c2
    else
        invariant = Universe(n)
    end

    #acceleration of the lead vehicle + time
    B = sparse([2], [1], [1.0], n, 1)
    U = Hyperrectangle(low=[-9.], high=[1.])
    c = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1.0]
    @system(x' = Ax + Bu + c, x ∈ invariant, u ∈ U)
end

# ## Hybrid System Setup

function platoon(; deterministic_switching::Bool=true,
                   c1=5.0,  # clock constraints
                   c2=5.0,  # clock constraints
                   tb=10.0,  # lower bound for loss of communication
                   tc=20.0, tr=20.0) # upper bound for loss of communication (tc) and reset time (tr)

    #three variables for each vehicle, (ei, d(et)/dt, ai) for
    #(spacing error, relative velocity, speed), and the last dimension is time
    n = 9 + 1
    var = @variables x₁, x₂, x₃, x₄, x₅, x₆, x₇, x₈, x₉, t

    #transition graph
    automaton = LightAutomaton(2)
    add_transition!(automaton, 1, 2, 1)
    add_transition!(automaton, 2, 1, 2)

    #modes
    mode1 = platoon_connected(deterministic_switching=deterministic_switching, c1=c1)
    mode2 = platoon_disconnected(deterministic_switching=deterministic_switching, c2=c2)
    modes = [mode1, mode2]

    #common reset
    reset = Dict(n => 0.)

    #transition l1 -> l2
    if deterministic_switching
        guard = Hyperplane(t == c1, var)
    else
        #tb <= t <= tc
        guard = HPolyhedron([tb <= t, t <= tc], var)
    end
    t1 = ConstrainedResetMap(n, guard, reset)

    #transition l2 -> l1
    if deterministic_switching
        guard = Hyperplane(t == c2, var)
    else
        guard = HalfSpace(t <= tr, var)
    end
    t2 = ConstrainedResetMap(n, guard, reset)
    resetmaps = [t1, t2]

    H = HybridSystem(automaton, modes, resetmaps, [AutonomousSwitching()])

    #initial condition is at the orgin in mode 1
    X0 = BallInf(zeros(n), 0.0)
    initial_condition = [(1, X0)]

    return IVP(H, initial_condition)
end

# ## Reachability settings
#
# The verification goal is to check whether the minimum distance between
# vehicles is preserved. The choice of the coordinate system is such that the
# minimum distance is a negative value.
# BNDxy Bounded time (no explicit bound on the number of transitions):
#    For all ``t ∈ [0, 20] [s]``,
#    ``x_1(t) ≥ −d_{min} [m], x_4(t) ≥ −d_{min} [m], x_7(t) ≥ −d_{min} [m]``,
#    where ``d_min = xy [m]``.
# BND50: ``d_{min} = 50``.
# BND42: ``d_{min} = 42``.
# BND30: ``d_{min} = 30``.
# UNBxy Unbounded time and unbounded switching: For all ``t ≥ 0 [s], x_1(t) ≥ −d_{min} [m]``,
# ``x_4(t) ≥ −d_{min} [m], x_7(t)``

# ## Results

# We use `LGG09` algorithm with ``δ=0.03`` and an octagonal template of
# directions.

octdirs = OctDirections(10);

#  ### PLAD01-BND30 (dense time)

prob_PLAD01_BND30 = platoon(; deterministic_switching=true);
imethod = TemplateHullIntersection(octdirs);
cmethod = LazyClustering(1);
alg = LGG09(δ=0.03, template=octdirs, approx_model=Forward(setops=octdirs));
sol_PLAD01_BND30 = solve(prob_PLAD01_BND30,
                         alg=alg,
                         clustering_method=cmethod,
                         intersection_method=imethod,
                         intersect_source_invariant=false,
                         tspan = (0.0 .. 20.0));

#--------------------------

using Plots

fig = plot();
plot!(fig, sol_PLAD01_BND30, vars=(0, 1), linecolor=:blue, color=:blue, alpha=0.8, lw=1.0,
    #tickfont=font(30, "Times"), guidefontsize=45, #!jl
    xlab="t",
    ylab="x1",
    xtick=[0, 5, 10, 15, 20.], ytick=[-30, -20, -10, 0],
    xlims=(0., 20.), ylims=(-31, 7))
fig
