using ReachabilityAnalysis, SparseArrays, JLD2

LazySets.set_ztol(Float64, 1e-14)

const x25 = [zeros(24); 1.0; zeros(23)]
const x25e = vcat(x25, 0.0)

examples_dir = normpath(@__DIR__, "..", "..", "..", "examples")
building_path = joinpath(examples_dir, "Building", "building.jld2")

function building_BLDF01()
    @load building_path H
    vars = collect(keys(H.ext[:variables]))
    A = state_matrix(mode(H, 1))
    n = size(A, 1) - 1
    A = A[1:n, 1:n]
    B = input_matrix(mode(H, 1))[1:n, 1]
    U = Hyperrectangle(low=[0.8], high=[1.0])
    S = @system(x' = Ax + Bu, u ∈ U, x ∈ Universe(n))

    #initial states
    center_X0 = [fill(0.000225, 10); fill(0.0, 38)]
    radius_X0 = [fill(0.000025, 10); fill(0.0, 14); 0.0001; fill(0.0, 23)]
    X0 = Hyperrectangle(center_X0, radius_X0)

    prob_BLDF01 = InitialValueProblem(S, X0)
end

function building_BLDC01()
    @load building_path H
    vars = collect(keys(H.ext[:variables]))
    A = state_matrix(mode(H, 1))
    n = size(A, 1) - 1
    A = A[1:n, 1:n]
    Ae = copy(transpose(hcat(transpose(hcat(A, zeros(48))), zeros(49))))
    S = LinearContinuousSystem(Ae)

    #initial states
    center_X0 = [fill(0.000225, 10); fill(0.0, 38)]
    radius_X0 = [fill(0.000025, 10); fill(0.0, 14); 0.0001; fill(0.0, 23)]
    X0 = Hyperrectangle(center_X0, radius_X0)
    U = Hyperrectangle(low=[0.8], high=[1.0])
    X0 = X0 * U

    prob_BLDC01 = InitialValueProblem(S, X0)
end

using Plots

prob_BLDF01 = building_BLDF01()

sol_BLDF01_dense = solve(prob_BLDF01, T=20.0, alg=LGG09(δ=0.004, template=[x25, -x25]));

plot(sol_BLDF01_dense, vars=(0, 25), linecolor=:blue, color=:blue, alpha=0.8, lw=1.0, xlab="t", ylab="x25")

ρ(x25, sol_BLDF01_dense)

ρ(x25, sol_BLDF01_dense) <= 5.1e-3 # BLDF01 - BDS01

ρ(x25, sol_BLDF01_dense) <= 4e-3 # BLDF01 - BDU01

ρ(x25, sol_BLDF01_dense(20.0))

ρ(x25, sol_BLDF01_dense(20.0)) <= -0.78e-3 # BLDF01 - BDU02

sol_BLDF01_discrete = solve(prob_BLDF01, T=20.0, alg=LGG09(δ=0.01, template=[x25, -x25], approx_model=NoBloating()));

plot(sol_BLDF01_discrete, vars=(0, 25), linecolor=:blue, color=:blue, alpha=0.8, lw=1.0, xlab="t", ylab="x25")

ρ(x25, sol_BLDF01_discrete)

ρ(x25, sol_BLDF01_discrete) <= 5.1e-3 # BLDF01 - BDS01

ρ(x25, sol_BLDF01_discrete)

ρ(x25, sol_BLDF01_discrete) <= 4e-3 # BLDF01 - BDU01

ρ(x25, sol_BLDF01_discrete(20.0))

ρ(x25, sol_BLDF01_discrete(20.0)) <= -0.78e-3 # BLDF01 - BDU02

prob_BLDC01 = building_BLDC01()

sol_BLDC01_dense = solve(prob_BLDC01, T=20.0, alg=LGG09(δ=0.006, template=[x25e, -x25e]))

plot(sol_BLDC01_dense, vars=(0, 25), linecolor=:blue, color=:blue, alpha=0.8, lw=1.0, xlab="t", ylab="x25")

ρ(x25e, sol_BLDC01_dense)

ρ(x25e, sol_BLDC01_dense) <= 5.1e-3 # BLDC01 - BDS01

ρ(x25e, sol_BLDC01_dense) <= 4e-3 # BLDC01 - BDU01

ρ(x25, sol_BLDF01_discrete(20.0))

ρ(x25e, sol_BLDC01_dense(20.0)) <= -0.78e-3 # BLDC01 - BDU02

sol_BLDC01_discrete = solve(prob_BLDC01, T=20.0, alg=LGG09(δ=0.01, template=[x25e, -x25e], approx_model=NoBloating()))

plot(sol_BLDC01_discrete, vars=(0, 25), linecolor=:blue, color=:blue, alpha=0.8, lw=1.0, xlab="t", ylab="x25")

ρ(x25e, sol_BLDC01_discrete)

ρ(x25e, sol_BLDC01_discrete) <= 5.1e-3 # BLDC01 - BDS01

ρ(x25e, sol_BLDC01_discrete)

ρ(x25e, sol_BLDC01_discrete) <= 4e-3 # BLDC01 - BDU01

ρ(x25e, sol_BLDC01_discrete(20.0))

ρ(x25e, sol_BLDC01_discrete(20.0)) <= -0.78e-3 # BLDC01 - BDU02

