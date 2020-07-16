# # Building
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/models/platoon.ipynb)
#
#md # !!! note "Overview"
#md #     System type: Affine system\
#md #     State dimension: 48\
#md #     Application domain: Building engineering
#
# ## Model description
#
# This benchmark is quite straightforward: The system is described by
# ``\dot{x}(t) = Ax(t) + Bu(t)``, ``u(t) \in \mathcal{U}``, ``y(t) = Cx(t)``
#
# Discrete-time analysis for the building system should use a step size of ``0.01``.
# There are two versions of this benchmark:
# - The inputs can change arbitrarily over time: $\forall t: u(t)\in \mathcal{U}$.
# - (constant inputs) The inputs are uncertain only in their initial value, and
#    constant over time: ``u(0)\in \mathcal{U}``, ``\dot u (t)= 0``. The purpose
#    of this model instance is to accommodate tools that cannot handle
#    time-varying inputs.

using ReachabilityAnalysis, SparseArrays, JLD2

LazySets.set_ztol(Float64, 1e-14)

const x25 = [zeros(24); 1.0; zeros(23)]
const x25e = vcat(x25, 0.0);
building_path = joinpath(@__DIR__, "building.jld2")

function building_BLDF01()
    @load building_path H
    vars = collect(keys(H.ext[:variables]))
    A = state_matrix(mode(H, 1))
    n = size(A, 1) - 1
    A = A[1:n, 1:n]
    B = input_matrix(mode(H, 1))[1:n, 1]
    U = Hyperrectangle(low=[0.8], high=[1.0])
    S = @system(x' = Ax + Bu, u ∈ U, x ∈ Universe(n));

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

# ## Reachability settings
#
# The verification goal is to check whether the displacement ``y_1`` of the top
# floor of the building remains below a given bound. In addition to the safety
# specification from the original benchmark, there are two UNSAT instances that
# serve as sanity checks to ensure that the model and the tool work as intended.
# But there is a caveat: In principle, verifying an UNSAT instance only makes
# sense formally if a witness is provided (counter-example, under-approximation,
# etc.). Since most of the participating tools do not have this capability, we
# run the tools with the same accuracy settings on an SAT-UNSAT pair of
# instances. The SAT instance demonstrates that the over-approximation is not
# too coarse, and the UNSAT instance indicates that the over-approximation is
# indeed conservative.
# - [BDS01] Bounded time, safe property: For all ``t \in [0, 20]``,
#    ``y_1(t) \leq 5.1\cdot 10^{-3}``. This property is assumed to be satisfied.
# - [BDU01] Bounded time, unsafe property: For all ``t \in [0, 20]``,
#    ``y_1(t) \leq 4\cdot 10^{-3}``. This property is assumed to be violated.
#    Property BDU01 serves as a sanity check. A tool should be run with the same
#    accuracy settings on BLDF01-BDS01 and BLDF01-BDU01, returning UNSAT on the
#    former and SAT on the latter.
# - [BDU02] Bounded time, unsafe property: The forbidden states are
#    ``\{ y_1(t) \leq -0.78\cdot 10^{-3} \wedge t = 20\}``. This property is
#    assumed to be violated for BLDF01 and satisfied for BLDC01. Property BDU02
#    serves as a sanity check to confirm that time-varying inputs are taken into
#    account. A tool should be run with the same accuracy settings on BLDF01-BDU02
#    and BLDC01-BDU02, returning UNSAT on the former and SAT on the latter.


# ## Results

using Plots

# ### BLDF01

prob_BLDF01 = building_BLDF01()

# ### Dense time
sol_BLDF01_dense = solve(prob_BLDF01, T=20.0, alg=LGG09(δ=0.004, template=x25))
plot(sol_BLDF01_dense, vars=(0, 1), linecolor=:blue, color=:blue, alpha=0.8, lw=1.0,
    xlab="t",
    ylab="x1")

# ### Discrete time
sol_BLDF01_discrete = solve(prob_BLDF01, T=20.0, alg=LGG09(δ=0.01, template=x25, approx_model=NoBloating()))
plot(sol_BLDF01_discrete, vars=(0, 1), linecolor=:blue, color=:blue, alpha=0.8, lw=1.0,
    xlab="t",
    ylab="x1")


# ### BLDC01

prob_BLDC01 = building_BLDC01()

# ### Dense time
sol_BLDC01_dense = solve(prob_BLDC01, T=20.0, alg=LGG09(δ=0.006, template=x25e))
plot(sol_BLDC01_dense, vars=(0, 1), linecolor=:blue, color=:blue, alpha=0.8, lw=1.0,
    xlab="t",
    ylab="x1")

# ### Discrete time
sol_BLDC01_discrete = solve(prob_BLDC01, T=20.0, alg=LGG09(δ=0.01, template=x25e, approx_model=NoBloating()))
plot(sol_BLDC01_discrete, vars=(0, 1), linecolor=:blue, color=:blue, alpha=0.8, lw=1.0,
    xlab="t",
    ylab="x1")
