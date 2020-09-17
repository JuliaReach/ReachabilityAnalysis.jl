# # Brusselator
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/models/Brusselator.ipynb)
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

# The numerical values for the model's constants (in their respective units) are
# given in the following table.
#
# |Quantity|Value|
# |-----|-------|
# |A    | 1     |
# |B    | 1.5   |

using ReachabilityAnalysis

const A = 1.0
const B = 1.5
const B1 = B + 1

@taylorize function brusselator!(du, u, p, t)
    x, y = u
    x² = x * x
    aux = x² * y
    du[1] = A + aux - B1*x
    du[2] = B*x - aux
    return du
end

#md # !!! tip "Performance tip"
#md #     The auxiliary variables `B1`, `x²` and `aux` have been defined to make
#md #     better use of `@taylorize` and help to reduce allocations.

# ## Reachability settings
#
# The initial set is defined by ``x \in [0.8, 1]``, ``y \in [0, 0.2]``.
# These settings are taken from [^CAS13].

U₀ = (0.8 .. 1.0) × (0.0 .. 0.2); #!jl
prob = @ivp(u' = brusselator!(U), u(0) ∈ U₀, dim: 2); #!jl

# ## Results

# We use `TMJets` algorithm with sixth-order expansion in time and second order expansion
# in the spatial variables.

sol = solve(prob, T=18.0, alg=TMJets(orderT=6, orderQ=2)); #!jl

#-

using Plots #!jl

solz = overapproximate(sol, Zonotope) #!jl
plot(solz, vars=(1, 2), xlab="x", ylab="y", lw=0.2, color=:blue, lab="Flowpipe", legend=:bottomright) #!jl
plot!(U₀, color=:orange, lab="Uo") #!jl

# We observe that the system converges to the equilibrium point `(1.0, 1.5)`.

# Below we plot the flowpipes projected into the time domain.

plot(solz, vars=(0, 1), xlab="t", lw=0.2, color=:blue, lab="x(t)", legend=:bottomright)
plot!(solz, vars=(0, 2), xlab="t", lw=0.2, color=:red, lab="y(t)")

# ## References

# [^CAS13]: X. Chen, E. Abraham, S. Sankaranarayanan. *Flow*: An Analyzer for Non-Linear Hybrid Systems.* In Proceedings of the 25th International Conference on Computer Aided Verification (CAV’13). Volume 8044 of LNCS, pages 258-263, Springer, 2013.
