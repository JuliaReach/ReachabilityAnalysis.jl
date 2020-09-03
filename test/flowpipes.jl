@testset "Flowpipe constructors" begin

    # uninitialized array constructor
    # if the set type is passed, then a ReachSet is used by default
    fp = Flowpipe(undef, Singleton{Float64, Vector{Float64}}, 3)
    @test length(fp) == 3
    @test_throws UndefRefError dim(fp) # fp is not empty but its dimension is not defined
    @test setrep(fp) == Singleton{Float64, Vector{Float64}}

    # the reach-set type can be passed directly
    fpr = Flowpipe(undef, ReachSet{Float64, Singleton{Float64, Vector{Float64}}}, 3)
    @test length(fpr) == 3
    @test typeof(fpr) == typeof(fp)
    @test setrep(fpr) == Singleton{Float64, Vector{Float64}}

    # try another reach-set type
    RSP = SparseReachSet{Float64, Singleton{Float64, Vector{Float64}}, 2}
    fp = Flowpipe(undef, RSP, 3)
    @test length(fp) == 3
    @test dim(fp) == 2 # in this set type, the dimension is known
    @test setrep(fp) == Singleton{Float64, Vector{Float64}}
end

@testset "Hybrid flowpipe time evaluation" begin
    X = Interval(0 .. 1)
    δ = 0.1
    F1 = Flowpipe([ReachSet(X, (0 .. δ) + k*δ) for k in 0:10])
    F2 = Flowpipe([ReachSet(X, (0 .. δ) + k*δ) for k in 9:20])
    F2′ = Flowpipe([ReachSet(X, (0 .. δ) + k*δ*0.9) for k in 9:20])
    H = HybridFlowpipe([F1, F2])
    H′ = HybridFlowpipe([F1, F2′])

    times = 0:0.01:1.9
    values = [[x] for x in (@. 1 - exp(-times))]
    @test all(vi ∈ H(ti) for (ti, vi) in zip(times, values))
    @test all(vi ∈ H′(ti) for (ti, vi) in zip(times, values))
end




using Reachability, Distributions, LazySets, ModelingToolkit

# ============================================================================
# Switched affine system
#
# Conider:
#     𝑥′(𝑡)=𝑓(𝑥(𝑡),𝑡)
#     x(𝑡0)=𝜙0
# where 𝑓(𝑥(𝑡),𝑡) is a periodic switched function, with period T = t₊ + t₋ such
# that 𝑓(𝑥(𝑡),𝑡) alternates between two modes:
# -> 𝑓(𝑥(𝑡),𝑡) = 𝑓₊(𝑥(𝑡),𝑡) = 𝛼₊*𝑥 + 𝛽₊, if t ∈ [t0, t0 + t₊]
# -> 𝑓(𝑥(𝑡),𝑡) = 𝑓₋(𝑥(𝑡),𝑡) = 𝛼₋*𝑥 + 𝛽₋, if t ∈ [t0 + t+, t0 + T]
#
# Here we consider the case of non-deterministic switching with jitter, i.e.
# uncertainty in the jumping time 𝛿, such that 𝑓 is equal to 𝑓₊
# during time t₊ ± 𝛿, and then changes to the mode 𝑓₋, where it stays
# during time t₋ ± 𝛿. After that the system jumpts again to the mode 𝑓₊
# and the cycle is repeated.

# The analytical solution is:
# -> 𝜙(𝑡) = (𝜙(𝑡0₊) + 𝛽₊/𝛼₊)*𝑒xp(𝛼₊*(𝑡 − 𝑡0₊)) + 𝛽₊/𝛼₊, if t ∈ [t0₊, t0₊ + t₊ ± 𝛿]
# -> 𝜙(𝑡) = (𝜙(𝑡0₋) + 𝛽₋/𝛼₋)*𝑒xp(𝛼₋*(𝑡 − 𝑡0₋)) + 𝛽₋/𝛼₋, if t ∈ [t0₋, t0₋ + t₋ ± 𝛿]
# where t0₊ is the time when the system swithched to 𝑓₊ mode and t0₋ is
# the time when the system jumped into the 𝑓₋ mode. The initial conditions at
# the begining of each mode 𝑓₊ and 𝑓₋ are 𝜙(𝑡0₊) and 𝜙(𝑡0₋) respectively.
# ============================================================================

@testset "Switched affine system" begin

    #We set the parameters 𝛼₊ = 𝛼₋ = -1.0, 𝛽₊ = 1.0 and 𝛽₋ = -𝛽₊ = -1.0:
    α₊ = -1.0
    β₊ = 1.0
    α₋ = -1.0
    β₋ = -1.0

    #The analytical solution for each mode:
    ϕ₊(t, t0, ϕt0) = exp(α₊*(t-t0))*(ϕt0 + β₊/α₊) - β₊/α₊
    ϕ₋(t, t0, ϕt0) = exp(α₋*(t-t0))*(ϕt0 + β₋/α₋) - β₋/α₋

    #Time limits of the mode f₊ (ξ₊ stands for the uncertainty in t₊, i.e. t₊+ξ₊)
    #In this example we choose ξ₊ ∈ [-𝛿, +𝛿] randomly asigned, since it's
    #non-deterministic. The case 𝛿 = 0 is the case of deterministic switching.
    left_times₊(k, t₊, t₋, ξ₋, ξ₊) = (k-1) * (t₊ + t₋) + sum(ξ₊[1:k-1] + ξ₋[1:k-1])
    right_times₊(k, t₊, t₋, ξ₋, ξ₊) = k*t₊ + (k-1)*t₋ + sum(ξ₊[1:k]) + sum(ξ₋[1:k-1])

    #Time limits of the mode f₋ (ξ₋ stands for the uncertainty in t₋, i.e. t₋+ξ₋)
    #Again, we will use ξ₋ ∈ [-𝛿, +𝛿]
    left_times₋(k, t₊, t₋, ξ₋, ξ₊) = right_times₊(k, t₊, t₋, ξ₋, ξ₊)
    right_times₋(k, t₊, t₋, ξ₋, ξ₊) = k* (t₊ + t₋) + sum(ξ₊[1:k] + ξ₋[1:k])

    #Analytical solution for N periods of the function f:
    function analytical_solution(N, ϕt0_ini, t₊, t₋, ξ₋, ξ₊; NPoints=3)

        #NPoints = number of sampling points in each mode
        @assert(NPoints>=2)

        ϕval = [[] for _ in 1:N]
        times = [[] for _ in 1:N]

        for k in 1:N

            # first semi-period (mode f₊)
            ϕt0 = k == 1 ? ϕt0_ini : ϕval[k-1][end]
            l = left_times₊(k, t₊, t₋, ξ₋, ξ₊)
            u = right_times₊(k, t₊, t₋, ξ₋, ξ₊)
            dt₊ = LinRange(l, u, NPoints) #dt₊ = range(l, u, step=Δ)
            t0 = l
            sol₊ = ϕ₊.(dt₊, t0, ϕt0)

            # second semi-period (mode f₋)
            l = left_times₋(k, t₊, t₋, ξ₋, ξ₊)
            u = right_times₋(k, t₊, t₋, ξ₋, ξ₊)
            dt₋ = LinRange(l, u, NPoints) #dt₋ = range(l, u, step=Δ)
            t0 = l
            ϕt0 = sol₊[end]
            sol₋ = ϕ₋.(dt₋, t0, ϕt0)

            times[k] = vcat(dt₊, dt₋)
            ϕval[k] = vcat(sol₊, sol₋)
        end
        return times, ϕval
    end

    function NRUNS_solutions(NRUNS, N, ϕt0_ini, t₊, t₋;
                             NPoints=3, δ=0.01, NPUNTOS=4)

        @assert NRUNS>1

        times = [[] for _  in 1:NRUNS]
        sol = [[] for _ in 1:NRUNS]

        I = LinRange(-δ, δ, NPUNTOS);
        t₊, t₋ = t₊, t₋

        #Solutions with random jumping times: ξ₋, ξ₊ ∈ [-𝛿, +𝛿]:
        for i in 1:NRUNS-2
            ξ₋, ξ₊ = rand(I, N), rand(I, N)
            times[i], sol[i] = analytical_solution(N, ϕt0_ini, t₊, t₋, ξ₋, ξ₊;
                                                   NPoints=NPoints)
        end

        #Extreme case 1: system always jump as soon as posible (ξ₋ = ξ₊ = -δ):
        i = NRUNS - 1
        I = LinRange(-δ, -δ, 1);
        ξ₋, ξ₊ = rand(I, N), rand(I, N)
        times[i], sol[i] = analytical_solution(N, ϕt0_ini, t₊, t₋, ξ₋, ξ₊;
                                               NPoints=NPoints)

        #Extreme case 2: system always jump late (ξ₋ = ξ₊ = +δ):
        i = NRUNS #i+=1
        I = LinRange(+δ, +δ, 1);
        ξ₋, ξ₊ = rand(I, N), rand(I, N)
        times[i], sol[i] = analytical_solution(N, ϕt0_ini, t₊, t₋, ξ₋, ξ₊;
                                               NPoints=NPoints)

        return times, sol
    end


    #Model of the Hybrid system to be solved with ReachabilityAnalysis:

    const vars = @variables x, a, b, t

    @taylorize function square_input₊!(du, u, p, t)
        x, α, β, t = u
        #Use: 𝛼 = 𝛼₊ and 𝛽 = 𝛽₊
        du[1] = α*x+β
        du[2] = zero(x)
        du[3] = zero(x)
        du[4] = one(t)
        return du
    end

    @taylorize function square_input₋!(du, u, p, t)
        x, α, β, t = u
        #Use: 𝛼 = 𝛼₊ and 𝛽 = 𝛽₊
        du[1] = α*x-β
        du[2] = zero(x)
        du[3] = zero(x)
        du[4] = one(t)
        return du
    end

    function hybrid_model(ϕt0_ini, α, β; deterministic_switching::Bool=true,
                          t₊=1.0, t₋=1.0, δ=0.01)

        n = 3 + 1;

        #variables:
        #x = 1;
        #𝛼 = 2;
        #𝛽 = 3;
        #t = 4;

        # initial set
        X0 = Hyperrectangle([ϕt0_ini, α, β, 0.0], [0.0, 0.0, 0.0, 0.0]);
        initial_condition = [(1, X0)]

        # discrete structure (graph)
        automaton = LightAutomaton(2) # 2 modes

        if deterministic_switching
            I₊ = LazySets.HalfSpace(t <= t₊, vars)
            I₋ = LazySets.HalfSpace(t <= t₋, vars)
        else
            I₊ = LazySets.HalfSpace(t <= t₊+δ, vars)#Universe(n)
            I₋ = LazySets.HalfSpace(t <= t₋+δ, vars)#Universe(n)
        end

        m1 = @system(x' = square_input₊!(x), dim: 4, x ∈ I₊);
        m2 = @system(x' = square_input₋!(x), dim: 4, x ∈ I₋);

        modes = [m1, m2]

        #common reset
        reset = Dict(n => 0.)

        # transition m1 -> m2
        add_transition!(automaton, 1, 2, 1)

        if deterministic_switching
            guard = Hyperplane(t == t₊, vars)
        else
            guard = HPolyhedron([t₊-δ <= t, t <= t₊+δ], vars) #TODO: t₊-δ <= t <= t₊+δ
        end
        t1 = ConstrainedResetMap(n, guard, reset)

        # transition m2 -> m1
        add_transition!(automaton, 2, 1, 2)

        if deterministic_switching
            guard = Hyperplane(t == t₋, vars)
        else
            guard = HPolyhedron([t₋-δ <= t, t <= t₋+δ], vars)
        end
        t2 = ConstrainedResetMap(n, guard, reset)

        # transition annotations
        resetmaps = [t1, t2];

        # switching
        switchings = [AutonomousSwitching()];

        H = HybridSystem(automaton, modes, resetmaps, switchings);

        return InitialValueProblem(H, initial_condition);

    end

    #Parameters and initial condicion for t=0:
    NRUNS = 3 #number of solutions: 1 random and the 2 extreme cases
    N = 10 # number of periods
    δ = 0.01;
    ϕt0_ini = 0.0;
    t₊, t₋ = 0.5, 1.0;
    boxdirs = BoxDirections(4)

    #Analytical solutions (NRUNS trajectories):
    times, sol = NRUNS_solutions(NRUNS, N, ϕt0_ini, t₊, t₋; δ=δ);

    #Flowpipe calculation from ReachabilityAnalysis:
    prob = hybrid_model(ϕt0_ini, α₊, β₊,
                        deterministic_switching=false, t₊=t₊, t₋=t₋, δ=δ);
                        tf = N*(t₊+t₋+2*δ);

    sol_ext = solve(prob,
        tspan=(0.0, N*(t₊+t₋+2*δ)),
        alg=TMJets(abs_tol=1e-9, orderT=4, orderQ=1, disjointness=BoxEnclosure()),
        #max_jumps=4,
        intersect_source_invariant=false,
        intersection_method=TemplateHullIntersection(boxdirs),#no necesita Polyhedra
        clustering_method=BoxClustering(1),#LazyClustering(1),
        disjointness_method=BoxEnclosure(),#ZonotopeEnclosure());
        fixpoint_check=false);

    #Verify if the analytical solution is inside the flowpipe:

    t2 = reduce(vcat,reduce(vcat, tiempos));
    v2 = reduce(vcat,reduce(vcat, sol));
    @test all(vi ∈ sol_ext(ti) for (ti, vi) in zip(t2, v2))

    sol_extz = overapproximate(sol_ext, Zonotope);
    @test all(vi ∈ sol_extz(ti) for (ti, vi) in zip(t2, v2))

end
