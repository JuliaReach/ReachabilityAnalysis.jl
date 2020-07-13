# # Platooning
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/models/Platoon.ipynb)
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

using ReachabilityAnalysis, SparseArrays

LazySets.set_ztol(Float64, 1e-15);

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

#-

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

#-

function platoon(; deterministic_switching::Bool=true,
                   c1=5.0,  # clock constraints
                   c2=5.0,  # clock constraints
                   tb=10.0,  # lower bound for loss of communication
                   tc=20.0, tr=20.0) # upper bound for loss of communication (tc) and reset time (tr)

    #three variables for each vehicle, (ei, d(et)/dt, ai) for
    #(spacing error, relative velocity, speed), and the last dimension is time
    n = 9 + 1

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
        guard = Hyperplane(sparsevec([n], [1.], n), c1) # t == c1
    else
        #tb <= t <= tc
        guard = HPolyhedron([HalfSpace(sparsevec([n], [-1.], n), -tb),
                             HalfSpace(sparsevec([n], [1.], n), tc)])
    end
    t1 = ConstrainedResetMap(n, guard, reset)

    #transition l2 -> l1
    if deterministic_switching
        guard = Hyperplane(sparsevec([n], [1.], n), c2) # t == c2
    else
        guard = HalfSpace(sparsevec([n], [1.], n), tr) # t <= tr
    end
    t2 = ConstrainedResetMap(n, guard, reset)
    resetmaps = [t1, t2]

    H = HybridSystem(automaton, modes, resetmaps, [AutonomousSwitching()])

    #initial condition is at the orgin in mode 1
    X0 = BallInf(zeros(n), 0.0)
    initial_condition = [(1, X0)]

    return IVP(H, initial_condition)
end

#=
function dmin_specification(sol, dmin)
    return (-ρ(sparsevec([1], [-1.0], 10), sol) > dmin) &&
           (-ρ(sparsevec([4], [-1.0], 10), sol) > dmin) &&
           (-ρ(sparsevec([7], [-1.0], 10), sol) > dmin)
end
=#

#-

# ## Results

# TODO: We use `LGG09` algorithm with sixth-order expansion in time and second order expansion
# in the spatial variables.

#boxdirs = BoxDirections{Float64, Vector{Float64}}(10); #!jl
octdirs = CustomDirections([Vector(vi) for vi in OctDirections(10)]);

# ----------------------------------------
#  PLAD01-BND30 (dense time)
# ----------------------------------------
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
#property = dmin_specification(sol_PLAD01_BND30, -30.0) #!jl

#-

using Plots, Plots.PlotMeasures, LaTeXStrings #!jl

fig = Plots.plot();
Plots.plot!(fig, sol_PLAD01_BND30, vars=(0, 1), linecolor=:blue, color=:blue, alpha=0.8, lw=1.0,
    #tickfont=font(30, "Times"), guidefontsize=45, #!jl
    xlab=L"t",
    ylab=L"x_{1}",
    xtick=[0, 5, 10, 15, 20.], ytick=[-30, -20, -10, 0],
    xlims=(0., 20.), ylims=(-31, 7),
    bottom_margin=6mm, left_margin=2mm, right_margin=4mm, top_margin=3mm)
    #,size=(1000, 1000)) #!jl
fig #!jl
