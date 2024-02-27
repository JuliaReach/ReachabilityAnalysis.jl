using ReachabilityAnalysis, JLD2, ReachabilityBase.CurrentPath
using ReachabilityBase.Arrays: SingleEntryVector

const x25 = SingleEntryVector(25, 48, 1.0)
const x25e = SingleEntryVector(25, 49, 1.0)

path = @current_path("Building", "building.jld2")

function building_common()
    @load path A B

    c = [fill(0.000225, 10); fill(0.0, 38)]
    r = [fill(0.000025, 10); fill(0.0, 14); 0.0001; fill(0.0, 23)]
    X0 = Hyperrectangle(c, r)

    U = Interval(0.8, 1.0)

    return A, B, X0, U
end

function building_BLDF01()
    A, B, X0, U = building_common()
    n = size(A, 1)
    S = @system(x' = A * x + B * u, u ∈ U, x ∈ Universe(n))
    prob_BLDF01 = InitialValueProblem(S, X0)
    return prob_BLDF01
end

using ReachabilityAnalysis: add_dimension

function building_BLDC01()
    A, B, X0, U = building_common()
    n = size(A, 1)
    Ae = add_dimension(A)
    Ae[1:n, end] = B
    prob_BLDC01 = @ivp(x' = Ae * x, x(0) ∈ X0 × U)
    return prob_BLDC01
end;

prob_BLDF01 = building_BLDF01();

sol_BLDF01_dense = solve(prob_BLDF01; T=20.0,
                         alg=LGG09(; δ=0.004, vars=(25), n=48))

fig = plot(sol_BLDF01_dense; vars=(0, 25), linecolor=:blue, color=:blue,
           alpha=0.8, lw=1.0, xlab="t", ylab="x25")

@assert ρ(x25, sol_BLDF01_dense) <= 5.1e-3 "the property should be proven"  # BLDF01 - BDS01
@assert !(ρ(x25, sol_BLDF01_dense) <= 4e-3) "the property should not be proven"  # BLDF01 - BDU01
ρ(x25, sol_BLDF01_dense)

@assert !(ρ(x25, sol_BLDF01_dense(20.0)) <= -0.78e-3) "the property should not be proven"  # BLDF01 - BDU02
ρ(x25, sol_BLDF01_dense(20.0))

sol_BLDF01_discrete = solve(prob_BLDF01; T=20.0,
                            alg=LGG09(; δ=0.01, vars=(25), n=48, approx_model=NoBloating()));

fig = plot(sol_BLDF01_discrete; vars=(0, 25), linecolor=:blue, color=:blue,
           alpha=0.8, lw=1.0, xlab="t", ylab="x25")

@assert ρ(x25, sol_BLDF01_discrete) <= 5.1e-3 "the property should be proven"  # BLDF01 - BDS01
@assert !(ρ(x25, sol_BLDF01_discrete) <= 4e-3) "the property should not be proven"  # BLDF01 - BDU01
ρ(x25, sol_BLDF01_discrete)

@assert !(ρ(x25, sol_BLDF01_discrete(20.0)) <= -0.78e-3) "the property should not be proven"  # BLDF01 - BDU02
ρ(x25, sol_BLDF01_discrete(20.0))

prob_BLDC01 = building_BLDC01();

sol_BLDC01_dense = solve(prob_BLDC01; T=20.0, alg=LGG09(; δ=0.005, vars=(25), n=49))

fig = plot(sol_BLDC01_dense; vars=(0, 25), linecolor=:blue, color=:blue,
           alpha=0.8, lw=1.0, xlab="t", ylab="x25")

@assert ρ(x25e, sol_BLDC01_dense) <= 5.1e-3 "the property should be proven"  # BLDC01 - BDS01
@assert !(ρ(x25e, sol_BLDC01_dense) <= 4e-3) "the property should not be proven"  # BLDC01 - BDU01
ρ(x25e, sol_BLDC01_dense)

@assert !(ρ(x25e, sol_BLDC01_dense(20.0)) <= -0.78e-3) "the property should not be proven"  # BLDC01 - BDU02
ρ(x25, sol_BLDF01_discrete(20.0))

sol_BLDC01_discrete = solve(prob_BLDC01; T=20.0,
                            alg=LGG09(; δ=0.01, vars=(25), n=49, approx_model=NoBloating()))

fig = plot(sol_BLDC01_discrete; vars=(0, 25), linecolor=:blue, color=:blue,
           alpha=0.8, lw=1.0, xlab="t", ylab="x25")

@assert ρ(x25e, sol_BLDC01_discrete) <= 5.1e-3 "the property should be proven" # BLDC01 - BDS01
@assert !(ρ(x25e, sol_BLDC01_discrete) <= 4e-3) "the property should not be proven" # BLDC01 - BDU01
ρ(x25e, sol_BLDC01_discrete)

@assert !(ρ(x25e, sol_BLDC01_discrete(20.0)) <= -0.78e-3) "the property should not be proven" # BLDC01 - BDU02
ρ(x25e, sol_BLDC01_discrete(20.0))
