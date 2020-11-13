```@meta
DocTestSetup = :(using ReachabilityAnalysis)
CurrentModule = ReachabilityAnalysis
```

# Benchmarks

## Scalability evaluation

Scalability is very important in the applicability of a tool. For illustration
purposes, in this section we consider the scalability of the hybrid reachability
solver, using different algorithm choices, for the filtered switched oscillator
model from [[FRE11]](@ref). The model consists of a two-dimensional switched oscillator
and a parametric number of filters which are used to *smooth* the oscilllator's state.
An interesting aspect of the model is that it is scalable: the total number of continuous
variables can be made arbitrarily large. Moreover, this is a challenging benchmark
since several dozens of reach-sets may take each discrete jumps, hence clustering methods are indispensible.

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
    a = LightAutomaton(4)
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


## Benchmark repository

The benchmark suite of JuliaReach is available at the
[ReachabilityBenchmarks](https://github.com/JuliaReach/ReachabilityBenchmarks) repository.

The repository includes:

- Benchmarks for all models in the `ReachabilityAnalysis.jl` test suite.
- Benchmarks for models which are not part of the `ReachabilityAnalysis.jl` test suite,
  and are used to detect regressions (i.e. examples that accidentally run slower,
  due to changes in core downstream libraries such as `LazySets.jl`).
- The [SLICOT models](http://slicot.org/matlab-toolboxes/model-reduction), from the
  SLICOT Model and Controller Reduction Toolbox, which reflect real world applications
  with dimensions ranging from several dozens to over 10.000.

## Repeatability evaluations

Traditionally, re-creation of computational results in research work is a challenging
task because details of the implementation are unavoidably absent in the paper.
Some authors post their code and data to their websites, but there is little formal
incentive to do so and no easy way to determine whether others can actually use the result.
As a consequence, computational results often become non reproducible -- even by the
research group which originally produced them -- after just a few years.

More recently, scientific conferences encourage and sometimes require that
authors improve the reproducibility of their computational results by providing suitable
*reproducibility evaluation* (RE) packages. Such RE are self-contained codes that
can be used to run and reproduce the results from the paper in a host machine.
Fortunately, the tools available for the Julia language are very convenenient
to prepare and distribute RE packages. The following RE packages are available:

- [ARCH2020 NLN RE](https://github.com/JuliaReach/ARCH2020_NLN_RE) -- Repeatability Evaluation package for the ARCH2020 Nonlinear Continuous and Hybrid Systems Competition.

- [ARCH2020 AFF RE](https://github.com/JuliaReach/ARCH2020_AFF_RE) -- Repeatability Evaluation package for the ARCH2020 Linear and Hybrid Systems Competition.

- [HSCC2019 RE](https://github.com/JuliaReach/HSCC2019_RE) -- Repeatability Evaluation (RE) package for the paper *JuliaReach: a Toolbox for Set-Based Reachability* published at the HSCC'2019 conference.

- [ARCH2019 RE](https://github.com/JuliaReach/ARCH2019_RE) -- Repeatability Evaluation package for the ARCH2019 Competition.

- [ARCH2018 RE](https://github.com/JuliaReach/ARCH2018_RE) -- Repeatability Evaluation package for the ARCH2018 Competition.

For installation instructions, see the `README` file in each package.
We have examples using virtual machines, and also using Docker containers for RE packages.
The advantage of using a Docker container is that downloading the necessary requirements
(including Julia itself) is an automated process.
