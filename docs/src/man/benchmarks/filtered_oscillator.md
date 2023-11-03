```@meta
DocTestSetup = :(using ReachabilityAnalysis)
CurrentModule = ReachabilityAnalysis
```

## Scalability evaluation

Scalability is very important in the applicability of a tool. For illustration
purposes, in this section we consider the scalability of the hybrid reachability
solver, using different algorithm choices, for the filtered switched oscillator
model from [[FRE11]](@ref). The model consists of a two-dimensional switched oscillator
and a parametric number of filters which are used to *smooth* the oscilllator's state.
An interesting aspect of the model is that it is scalable: the total number of continuous
variables can be made arbitrarily large. Moreover, this is a challenging benchmark
since several dozens of reach-sets may take each discrete jumps, hence clustering methods are indispensable.

The continuous variables ``x`` and ``y`` are used to denote the state of the oscillator,
and the remaining ``m`` variables are used for the state of the filters.
In this evaluation we consider that the number of filters ranges from 64 to 1024.
To measure the quality of the approximations, we consider the safety property given
by ``y(t) < 0.5`` for all ``t ∈ [0, T]``.

The filtered oscillator was also studied in [[BFFPS19]](@ref) to test the scalability
of a new scheme that exploits the sparsity of the hybrid automaton to efficiently
compute flowpipe-guard intersections. Such scheme is not considered in this section.

Recall that for a system of $n$ variables, a box overapproximation requires to
compute $n$ samples of the approximation boundary, whereas an octagon requires
$O(n^2)$ many. Hence, using an octagon template, while gives tighter results
in general, may incur in much higher computational cost if $n$ is high.

### Model

```@example filtered_oscillator
using ReachabilityAnalysis

one_loop_iteration = false
n0 = 4
n1 = (one_loop_iteration ? n0 + 1 : n0)
n = n1 + 2
z = zeros(n1)

## common flow
A = zeros(n, n)
A[1,1], A[2,2] = -2., -1.
A[3,1], A[3,3] = 5., -5.
for i = 4 : n-1
    A[i,i-1], A[i,i] = 5., -5.
end

function mode1(z)
    b = [1.4; -0.7; z]
    X = HPolyhedron([HalfSpace([-0.714286; -1.0; z], 0.0),  # 0.714286*x + y >= 0
                     HalfSpace([1.0; 0.0; z], 0.0)])  # x <= 0
    @system(x' = Ax + b, x ∈ X)
end

function mode2(z)
    b = [-1.4; 0.7; z]
    X = HPolyhedron([HalfSpace([1.0; 0.0; z], 0.0),  # x <= 0
                     HalfSpace([0.714286; 1.0; z], 0.0)])  # 0.714286*x + y <= 0
    @system(x' = Ax + b, x ∈ X)
end

function mode3(z)
    b = [1.4; -0.7; z]
    X = HPolyhedron([HalfSpace([-1.0; 0.0; z], 0.0),  # x >= 0
                     HalfSpace([-0.714286; -1.0; z], 0.0)])  # 0.714286*x + y >= 0
    @system(x' = Ax + b, x ∈ X)
end

function mode4(z, one_loop_iteration)
    b = [-1.4; 0.7; z]
    X = HPolyhedron([HalfSpace([0.714286; 1.0; z], 0.0),  # 0.714286*x + y <= 0
                     HalfSpace([-1.0; 0.0; z], 0.0)])  # x >= 0

    if one_loop_iteration
        ## k <= 2 (2.1 to avoid numerical issues)
        addconstraint!(X, HalfSpace([zeros(n-1); 1.], 2.1))
    end
    @system(x' = Ax + b, x ∈ X)
end


function filtered_oscillator_hybrid(n0, one_loop_iteration)

    n1 = (one_loop_iteration ? n0 + 1 : n0)
    n = n1 + 2
    z = zeros(n1)

    ## transition graph (automaton)
    a = GraphAutomaton(4)
    add_transition!(a, 3, 4, 1)
    add_transition!(a, 4, 2, 2)
    add_transition!(a, 2, 1, 3)
    add_transition!(a, 1, 3, 4)

    mode1 = mode1(z)
    mode2 = mode2(z)
    mode3 = mode3(z)
    mode4 = mode4(z, one_loop_iteration)
    m = [mode1, mode2, mode3, mode4]

    ## transitions

    ## transition l3 -> l4
    X_l3l4 = HPolyhedron([HalfSpace([-1.0; 0.0; z], 0.0),  # x >= 0
                          HalfSpace([-0.714286; -1.0; z], 0.0),  # 0.714286*x + y >= 0
                          HalfSpace([0.714286; 1.0; z], 0.0)])  # 0.714286*x + y <= 0
    if one_loop_iteration
        A_trans_34 = Matrix(1.0I, n, n)
        A_trans_34[n, n] = 2.  # k' = k * 2
        r1 = ConstrainedLinearMap(A_trans_34, X_l3l4)
    else
        r1 = ConstrainedIdentityMap(n, X_l3l4)
    end

    ## transition l4 -> l2
    X_l4l2 = HPolyhedron([HalfSpace([0.714286; 1.0; z], 0.0),  # 0.714286*x + y <= 0
                          HalfSpace([-1.0; 0.0; z], 0.0),  # x >= 0
                          HalfSpace([1.0; 0.0; z], 0.0)])  # x <= 0
    r2 = ConstrainedIdentityMap(n, X_l4l2)

    ## transition l2 -> l1
    X_l2l1 = HPolyhedron([HalfSpace([1.0; 0.0; z], 0.0),  # x <= 0
                          HalfSpace([-0.714286; -1.0; z], 0.0),  # 0.714286*x + y >= 0
                          HalfSpace([0.714286; 1.0; z], 0.0)])  # 0.714286*x + y <= 0
    r3 = ConstrainedIdentityMap(n, X_l2l1)

    ## transition l1 -> l3
    X_l1l3 = HPolyhedron([HalfSpace([-0.714286; -1.0; z], 0.0),  # 0.714286*x + y >= 0
                          HalfSpace([-1.0; 0.0; z], 0.0),  # x >= 0
                          HalfSpace([1.0; 0.0; z], 0.0)])  # x <= 0
    r4 = ConstrainedIdentityMap(n, X_l1l3)

    r = [r1, r2, r3, r4]

    ## switchings
    s = [HybridSystems.AutonomousSwitching()]

    return HybridSystem(a, m, r, s)
end

function fosc(X0; n0::Int=4, one_loop_iteration::Bool=false)
    H = filtered_oscillator_hybrid(n0, one_loop_iteration)
    return IVP(H, [(1, X0)])
end
```

## Results
