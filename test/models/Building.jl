using ReachabilityAnalysis, SparseArrays, JLD2

LazySets.set_ztol(Float64, 1e-14)

const x25 = [zeros(24); 1.0; zeros(23)]
const x25e = vcat(x25, 0.0)

examples_dir = normpath(@__DIR__, "..", "..", "..", "examples")
building_path = joinpath(examples_dir, "Building", "building.jld2")

function building_BLDF01()
    @load building_path A B
    n = size(A, 1)
    U = Interval(0.8, 1.0)
    S = @system(x' = Ax + Bu, u ∈ U, x ∈ Universe(n))

    # initial states
    center_X0 = [fill(0.000225, 10); fill(0.0, 38)]
    radius_X0 = [fill(0.000025, 10); fill(0.0, 14); 0.0001; fill(0.0, 23)]
    X0 = Hyperrectangle(center_X0, radius_X0)

    prob_BLDF01 = InitialValueProblem(S, X0)
end

using ReachabilityAnalysis: add_dimension

function building_BLDC01()
    @load building_path A B
    n = size(A, 1)
    U = Interval(0.8, 1.0)

    # initial states
    center_X0 = [fill(0.000225, 10); fill(0.0, 38)]
    radius_X0 = [fill(0.000025, 10); fill(0.0, 14); 0.0001; fill(0.0, 23)]
    X0 = Hyperrectangle(center_X0, radius_X0)

    Ae = add_dimension(A)
    Ae[1:n, end] = B
    prob_BLDC01 = @ivp(x' = Ae * x, x(0) ∈ X0 × U)
end

using Plots

prob_BLDF01 = building_BLDF01()

sol_BLDF01_dense = solve(prob_BLDF01, T=20.0, alg=LGG09(δ=0.004, vars=(25), n=48));

#plot(sol_BLDF01_dense, vars=(0, 25), linecolor=:blue, color=:blue, alpha=0.8, lw=1.0, xlab="t", ylab="x25")

ρ(x25, sol_BLDF01_dense)

ρ(x25, sol_BLDF01_dense) <= 5.1e-3 # BLDF01 - BDS01

ρ(x25, sol_BLDF01_dense) <= 4e-3 # BLDF01 - BDU01

ρ(x25, sol_BLDF01_dense(20.0))

ρ(x25, sol_BLDF01_dense(20.0)) <= -0.78e-3 # BLDF01 - BDU02

sol_BLDF01_discrete = solve(prob_BLDF01, T=20.0, alg=LGG09(δ=0.01, vars=(25), n=48, approx_model=NoBloating()));

#plot(sol_BLDF01_discrete, vars=(0, 25), linecolor=:blue, color=:blue, alpha=0.8, lw=1.0, xlab="t", ylab="x25")

ρ(x25, sol_BLDF01_discrete)

ρ(x25, sol_BLDF01_discrete) <= 5.1e-3 # BLDF01 - BDS01

ρ(x25, sol_BLDF01_discrete)

ρ(x25, sol_BLDF01_discrete) <= 4e-3 # BLDF01 - BDU01

ρ(x25, sol_BLDF01_discrete(20.0))

ρ(x25, sol_BLDF01_discrete(20.0)) <= -0.78e-3 # BLDF01 - BDU02

prob_BLDC01 = building_BLDC01()

sol_BLDC01_dense = solve(prob_BLDC01, T=20.0, alg=LGG09(δ=0.005, vars=(25), n=49))

#plot(sol_BLDC01_dense, vars=(0, 25), linecolor=:blue, color=:blue, alpha=0.8, lw=1.0, xlab="t", ylab="x25")

ρ(x25e, sol_BLDC01_dense)

ρ(x25e, sol_BLDC01_dense) <= 5.1e-3 # BLDC01 - BDS01

ρ(x25e, sol_BLDC01_dense) <= 4e-3 # BLDC01 - BDU01

ρ(x25, sol_BLDF01_discrete(20.0))

ρ(x25e, sol_BLDC01_dense(20.0)) <= -0.78e-3 # BLDC01 - BDU02

sol_BLDC01_discrete = solve(prob_BLDC01, T=20.0, alg=LGG09(δ=0.01, vars=(25), n=49, approx_model=NoBloating()))

#plot(sol_BLDC01_discrete, vars=(0, 25), linecolor=:blue, color=:blue, alpha=0.8, lw=1.0, xlab="t", ylab="x25")

ρ(x25e, sol_BLDC01_discrete)

ρ(x25e, sol_BLDC01_discrete) <= 5.1e-3 # BLDC01 - BDS01

ρ(x25e, sol_BLDC01_discrete)

ρ(x25e, sol_BLDC01_discrete) <= 4e-3 # BLDC01 - BDU01

ρ(x25e, sol_BLDC01_discrete(20.0))

ρ(x25e, sol_BLDC01_discrete(20.0)) <= -0.78e-3 # BLDC01 - BDU02

