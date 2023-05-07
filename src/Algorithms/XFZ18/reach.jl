using MultivariatePolynomials,
      DynamicPolynomials,
      JuMP,
      PolyJuMP,
      SumOfSquares,
      MathOptInterfaceMosek,
      SemialgebraicSets
using Reachability
using MathematicalSystems: PolynomialContinuousSystem, InitialValueProblem
const ∂ = differentiate

export XFZ18

"""
    XFZ18 <: ContinuousPost

Implementation of the reachability algorithm for the class of polynomial ODEs
with uncertain initial states (see [1], abbreviated `XFZ18`).
This method consists of reducing the Hamilton-Jacobi-Bellman equation to a
hierarchy of semidefinite programs that are solved using an SDP solver.

### Fields

- `options` -- an `Options` structure that holds the algorithm-specific options

### Algorithm

We refer to [1] for technical details.

- [1] *Xue, B., Fränzle, M., & Zhan, N. (2018, April). Under-Approximating Reach
Sets for Polynomial Continuous Systems. In [Proceedings of the 21st International
Conference on Hybrid Systems: Computation and Control (part of CPS Week)
(pp. 51-60). ACM.](https://dl.acm.org/citation.cfm?id=3178133)*
"""
struct XFZ18 <: ContinuousPost
    options::Options

    function XFZ18(𝑂::Options)

        # merge!(𝑂, :relaxation_degree)

        return new(𝑂)
    end
end

# build the SOS problem
function build_sos(𝒮, opt)
    T = opt[:T]

    # scale dynamics
    T = opt[:T]
    f = T * 𝒮.s.p

    # define polynomial symbolic variables
    vars = variables(𝒮.s) # TODO: variables(𝒮)
    @polyvar t

    k = opt[:relaxation_degree]

    # monomial vector up to order k
    # 0 <= sum_i alpha_i <= k, if alpha_i is the exponent of x_i
    X = monomials(vars, 0:k)
    XT = monomials([vars; t], 0:k)

    # create a SOS JuMP model to solve with Mosek
    solver = opt[:solver]
    model = SOSModel(with_optimizer(solver))

    # add unknown Φ to the model
    @variable(model, Φ, Poly(XT))

    # jacobian
    ∂t = α -> ∂(α, t)
    ∂xf = α -> ∂(α, x₁) * f[1] + ∂(α, x₂) * f[2]
    LΦ = ∂t(Φ) + ∂xf(Φ)

    # Φ(x, t) at time 0
    Φ₀ = subs(Φ, t => 0.0)

    # scalar variable
    @variable(model, ϵ)

    dom1 = @set t * (T - t) >= 0 && g >= 0
    dom2 = @set g >= 0
    @constraint(model, ϵ >= 0.0)
    @constraint(model, LΦ ∈ SOSCone(), domain = dom1)
    @constraint(model, ϵ - LΦ ∈ SOSCone(), domain = dom1)
    @constraint(model, Φ₀ - V₀ ∈ SOSCone(), domain = dom2)
    @constraint(model, ϵ + V₀ - Φ₀ ∈ SOSCone(), domain = dom2)

    @objective(model, Min, ϵ)
    return model
end

# solve model, check feasibility and return polynomials
function solve_sos(model; verbose=true)
    optimize!(model)

    if verbose
        println("Relaxation order : k = $k")
        println("JuMP.termination_status(model) = ", JuMP.termination_status(model))
        println("JuMP.primal_status(model) = ", JuMP.primal_status(model))
        println("JuMP.dual_status(model) = ", JuMP.dual_status(model))
        println("JuMP.objective_bound(model) = ", JuMP.objective_bound(model))
        println("JuMP.objective_value(model) = ", JuMP.objective_value(model))
    end

    # TODO: error if it fails?
end

function extract_approximations(model, 𝑂)

    # time horizon TODO : check consistency w/rescaling
    T = 𝑂[:T]

    # Recovering the solution:
    ϵopt = JuMP.objective_value(model)

    # Punder <= 0  TODO: @set Punder <= 0 ?
    Punder = subs(JuMP.value(model[:Φ]), t => T)

    # Pover <= 0   TODO: @set Pover <= 0 ?
    Pover = subs(JuMP.value(model[:Φ]), t => T) - ϵopt * (T + 1)

    return (ϵopt, Punder, Pover)
end

function post(𝒫::XFZ18, 𝒮::AbstractSystem, 𝑂::Options)

    # dynamics
    @assert 𝒮.s isa PolynomialContinuousSystem

    # construct sum-of-squares problem
    model = build_sos(𝒫, 𝒮, 𝑂)

    # solve the sum-of-squares optimization
    solve_sos(model, 𝑂)

    # extract under and over approximations
    (ϵopt, Punder, Pover) = extract_approximations(model)

    # returns the polynomial under and overapproximations of the reach set
    # for any t ∈ [0, T]
    return (Punder, Pover)
end
