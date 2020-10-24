# # Production-Destruction
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/models/ProductionDestruction.ipynb)
#
#md # !!! note "Overview"
#md #     System type: polynomial continuous system\
#md #     State dimension: 3\
#md #     Application domain: Chemical kinetics
#
# ## Model description
#
# A production-destruction system consists of an ordinary differential equation and two constraints:
# positivity and conservativity. It means that the system states are quantities which are always positive
# and the sum of these quantities is constant. This particular family of systems is often used to test
# the stability of integration schemes.
#
# As proposed by Michaelis-Menten theory [^KM18], a model for three quantities can be defined as follows:
#
# ```math
#    \begin{array}{lcl}
#    \dot{x} &=& \frac{-xy}{1+x}\\
#    \dot{y} &=& \frac{xy}{1+x} - a y\\
#    \dot{z} &=& ay
#    \end{array}
# ```
#
# with ``a=0.3`` and the initial condition ```x(0) = 9.98``, ``y(0) = 0.01`` and
# ``z(0) = 0.01``. In this model, ``x`` is the nutrients, ``y`` the phytoplankton
# and ``z`` the detritus.
#
# The constraints are:
#
# - ``x(t)``, ``y(t)``, ``z(t)`` are positive, and
# - ``x(t)+y(t)+z(t)=10`` for all ``t``.

# ## Analysis
#
# We are interested in computing the reachable tube till ``t=100``, and to verify
# (or embed) the  constraints. Three setups are considered, depending on the source
# of bounded uncertainties which is taken into account:
#
# - ``I``: ``x(0) \in [9.5, 10.0]``, i.e., uncertainty on the initial condition;
# - ``P``: ``a \in [0.296, 0.304]``, i.e., uncertainty on the parameter;
# - ``I \& P``: ``x(0) \in [9.7, 10.0]`` and ``a \in [0.298, 0.302]``, i.e., both uncertainties are mixed.
#
# In terms of objectives, at ``t=100``, the constraints ``10 \in x+y+z`` and
# ``x,y,z \geq 0`` have to be verified.

# Two measures are provided for comparison: the volume of the box (``x \times y \times z``)
# enclosing the final state (at ``t = 100``) and the total time of computation for
# evolution and verification. All of these results are obtained for the three setups.

# ## Results

# We observe that the specification is satisfied.
#
# Plots of ``z`` (in the ``[0, 11]`` range) w.r.t. time (in the ``[0, 100]`` range)
# for the three cases are overlaid in a single figure.

# there is a wide degree of difference in the volume of the final enclosure.
# However from Fig.~\ref{fig:productiondestruction} it is apparent that, except for the very tight result in DynIbex,
# similar errors are obtained for $z$ hence the differences in volume are largely due to convergence quality with respect
# to $x$,$y$, i.e., towards zero. Given these results, improved metrics would distinguish $z$ from $x$ and $y$.

using ReachabilityAnalysis, ModelingToolkit
using Plots, Plots.PlotMeasures, LaTeXStrings

@variables x y z
const positive_orthant = HPolyhedron([x >= 0, y >= 0, z >= 0], [x, y, z])

#md # !!! tip "Performance tip"
#md #     The functions below defines the system using some auxiliary variables to
#md #     get the most out of the `@taylorize` macro in terms of reducing allocations.
#md #     That said, defining `du[1] = -x*y / (1+x)` is also fine, but probably slower.

@taylorize function prod_dest_1!(du, u, params, t)
    local a = 0.3
    x, y, z = u[1], u[2], u[3]

    num = x * y
    den = 1 + x
    aux = num/den
    aux2 = a * y
    du[1] = -aux
    du[2] = aux - aux2
    du[3] = aux2
    return du
end

@taylorize function prod_dest_2!(du, u, params, t)
    x, y, z, a = u[1], u[2], u[3], u[4]

    num = x * y
    den = 1 + x
    aux = num/den
    aux2 = a * y
    du[1] = -aux
    du[2] = aux - aux2
    du[3] = aux2
    du[4] = zero(x)
    return du
end

function production_destruction(; case="I")
    if case == "I"
        X0 = (9.5 .. 10.0) × (0.01 .. 0.01) × (0.01 .. 0.01)
        prob = @ivp(x'= prod_dest_1!(x), dim:3, x(0) ∈ X0)

    elseif case == "P"
        X0 = (9.98 .. 9.98) × (0.01 .. 0.01) × (0.01 .. 0.01) × (0.296 .. 0.304)
        prob = @ivp(x'= prod_dest_2!(x), dim:4, x(0) ∈ X0)

    elseif case == "IP"
        X0 = (9.7 .. 10.0) × (0.01 .. 0.01) × (0.01 .. 0.01) × (0.298 .. 0.302)
        prob = @ivp(x'= prod_dest_2!(x), dim:4, x(0) ∈ X0)

    else
        error("case = $case is not implemented")
    end
    return prob
end

function prod_dest_property(solz)
    X = project(solz(100.0), vars=(1, 2, 3))

    ## check that all variables are nonnegative
    nonnegative = X ⊆ positive_orthant

    ## compute the volume of the last reach-set
    H = overapproximate(X, Hyperrectangle)
    vol = volume(H)

    ## check that that 10.0 belongs to the minkowski sum of the reach-sets projected in each coordinate
    B = convert(IntervalBox, H) # get the product-of-intervals representation
    contains_10 = 10 ∈ sum(B)

    return nonnegative && contains_10, vol
end

# For all three cases we use ``n_T = 7``, ``n_Q = 1``, and an adaptive absolute
# tolerance with initial value ``10^{-11}`` (resp ``10^{-12}``) is used for ``I``
# and ``I\&P`` (resp. ``P``).

# ### Case ``I``: uncertain initial states

prob = production_destruction(case="I")
alg = TMJets(abs_tol=1e-11, orderT=7, orderQ=1, adaptive=true)
sol_pd1 = solve(prob, T=100.0, alg=alg)
sol_pd1z = overapproximate(sol_pd1, Zonotope);

# Verifying that the specification holds:
property, vol = prod_dest_property(sol_pd1z)
validation = Int(property)
final_volume = trunc(vol, sigdigits=2)
println("volume final box, case P : $final_volume")

# ### Case ``P``: uncertain parameter

prob = production_destruction(case="P")
alg = TMJets(abs_tol=1e-12, orderT=7, orderQ=1, adaptive=true)
sol_pd2 = solve(prob, T=100.0, alg=alg)
sol_pd2z = overapproximate(sol_pd2, Zonotope);

# Verifying that the specification holds:
property, vol = prod_dest_property(sol_pd2z)
validation = Int(property)
final_volume = trunc(vol, sigdigits=2)
println("volume final box, case P : $final_volume")

# ### Case ``I \& P``: uncertain initial states and parameters

prob = production_destruction(case="IP")
alg = TMJets(abs_tol=1e-11, orderT=7, orderQ=1, adaptive=true)
sol_pd3 = solve(prob, T=100.0, alg=alg)
sol_pd3z = overapproximate(sol_pd3, Zonotope);

# verify that specification holds
property, vol = prod_dest_property(sol_pd3z)
validation = Int(property)
final_volume = trunc(vol, sigdigits=2)
println("volume final box, case I&P : $final_volume")

#-

# Plots of ``z`` (in the ``[0, 11]`` range) w.r.t. time (in the ``[0, 100]`` range)
# for the three cases are overlaid in a single figure.

fig = Plots.plot()
dt = 0 .. 100

Plots.plot!(fig, sol_pd3z(dt),  vars=(0, 3), linecolor="red", color=:red, alpha=3.0, leg=:bottomright, lab="I & P")

Plots.plot!(fig, sol_pd2z(dt), vars=(0, 3), linecolor="blue", color=:blue, alpha=0.8,
    tickfont=font(30, "Times"), guidefontsize=45,
    xlab=L"t",
    ylab=L"z",
    xtick=[0., 25., 50., 75., 100.], ytick=[0.0, 2.5, 5.0, 7.5, 10.0],
    xlims=(0., 100.5), ylims=(0.0, 11.0),
    bottom_margin=6mm, left_margin=2mm, right_margin=6mm, top_margin=3mm,
    size=(1000, 1000), lab="P", legendfontsize=20)

Plots.plot!(fig, sol_pd1z(dt), vars=(0, 3), linecolor="yellow", color=:yellow, alpha=0.3, lab="I")

# ## References

# [^KM18]: Kopecz, Stefan, and Andreas Meister. *On order conditions for modified Patankar–Runge–Kutta schemes.* Applied Numerical Mathematics 123 (2018): 159-179.

# [^CAS13]: X. Chen, E. Abraham, S. Sankaranarayanan. *Flow*: An Analyzer for Non-Linear Hybrid Systems.* In Proceedings of the 25th International Conference on Computer Aided Verification (CAV’13), Volume 8044 of LNCS, pages 258-263, Springer, 2013.
