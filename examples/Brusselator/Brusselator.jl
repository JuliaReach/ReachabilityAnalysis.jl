# # Brusselator
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

using ReachabilityAnalysis #!jl

const A = 1.0
const B = 1.5
const B1 = B + 1

@taylorize function brusselator!(du, u, p, t)
    x, y = u
    x² = x * x
    aux = x² * y
    du[1] = A + aux - B1 * x
    du[2] = B * x - aux
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
prob = @ivp(u' = brusselator!(U), u(0) ∈ U₀, dim:2); #!jl

# ## Results

# We use `TMJets` algorithm with sixth-order expansion in time and second order expansion
# in the spatial variables.

sol = solve(prob; T=18.0, alg=TMJets20(; orderT=6, orderQ=2)); #!jl

#-

using Plots #!jl

fig = plot(sol; vars=(1, 2), xlab="x", ylab="y", lw=0.2, color=:blue, lab="Flowpipe", #!jl
           legend=:bottomright) #!jl
plot!(U₀; color=:orange, lab="Uo") #!jl

#!jl import DisplayAs  #hide
#!jl DisplayAs.Text(DisplayAs.PNG(fig))  #hide

# We observe that the system converges to the equilibrium point `(1.0, 1.5)`.

# Below we plot the flowpipes projected into the time domain.

fig = plot(sol; vars=(0, 1), xlab="t", lw=0.2, color=:blue, lab="x(t)", legend=:bottomright)  #!jl
plot!(sol; vars=(0, 2), xlab="t", lw=0.2, color=:red, lab="y(t)")  #!jl

#!jl DisplayAs.Text(DisplayAs.PNG(fig))  #hide

# ## Changing the initial volume

# The model was considered in [^GCLASG20] but using a different set of initial conditions. Let us parametrize
# the initial states as a ball centered at ``x = y = 1`` and radius ``r > 0``:

U0(r) = Singleton([1.0, 1.0]) ⊕ BallInf(zeros(2), r)

#-

# The parametric initial-value problem is defined accordingly.

bruss(r) = @ivp(u' = brusselator!(u), u(0) ∈ U0(r), dim:2)

# First we solve for ``r = 0.01``:

sol_01 = solve(bruss(0.01); T=30.0, alg=TMJets20(; orderT=6, orderQ=2))  #!jl

LazySets.set_ztol(Float64, 1e-15)  #!jl

fig = plot(sol_01; vars=(1, 2), xlab="x", ylab="y", lw=0.2, color=:blue, lab="Flowpipe (r = 0.01)",  #!jl
           legend=:bottomright)  #!jl

plot!(U0(0.01); color=:orange, lab="Uo", xlims=(0.6, 1.3))  #!jl

#!jl DisplayAs.Text(DisplayAs.PNG(fig))  #hide

# We observe that the wrapping effect is controlled and the flowpipe doesn't blow up even for the large time horizon ``T = 30.0``.
# Next we plot the flowpipe zoomed to the last portion and compare ``r = 0.01`` with a set of larger initial states, ``r = 0.1``.

sol_1 = solve(bruss(0.1); T=30.0, alg=TMJets20(; orderT=6, orderQ=2)) #!jl

fig = plot(; xlab="x", ylab="y", xlims=(0.9, 1.05), ylims=(1.43, 1.57), legend=:bottomright)  #!jl

plot!(sol_1; vars=(1, 2), lw=0.2, color=:red, lab="r = 0.1", alpha=0.4)  #!jl

plot!(sol_01; vars=(1, 2), lw=0.2, color=:blue, lab="r = 0.01")  #!jl

#!jl DisplayAs.Text(DisplayAs.PNG(fig))  #hide

# The volume at time ``T = 9.0`` can be obtained by evaluating the flowpipe and computing the volume of the hyperrectangular overapproximation:

vol_01 = volume(set(overapproximate(sol_01(9.0), Hyperrectangle)))  #!jl

#-

vol_1 = volume(set(overapproximate(sol_1(9.0), Hyperrectangle)))  #!jl

# ## References

# [^CAS13]: X. Chen, E. Abraham, S. Sankaranarayanan. *Flow*: An Analyzer for Non-Linear Hybrid Systems.* In Proceedings of the 25th International Conference on Computer Aided Verification (CAV’13). Volume 8044 of LNCS, pages 258-263, Springer, 2013.

# [^GCLASG20]: S. Gruenbacher, J. Cyranka, M. Lechner, Md. Ariful Islam,, Scott A. Smolka, R. Grosu *Lagrangian Reachtubes: The Next Generation*. arXiv: 2012.07458
