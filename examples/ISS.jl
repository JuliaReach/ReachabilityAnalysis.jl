# # International Space Station
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/models/ISS.ipynb)
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

using ReachabilityAnalysis, JLD2
using ReachabilityAnalysis: add_dimension

LazySets.set_ztol(Float64, 1e-15); #!jl
ISS_path = joinpath(@__DIR__, "ISS.jld2"); #!jl

@load ISS_path C; #!jl
const C3 = C[3, :]; #!jl # variable y₃
const C3_ext = vcat(C3, fill(0.0, 3)); #!jl

# ## Reachability settings
#
# TODO:
# ISSC01 description
# The initial set ``X_0`` is defined by ``B = {x_i \in [-0001,+0.0001]``, ``\forall i}`` and
# ``U \in [0, 0.1] \times [0.8, 1.] \times [0.9, 1.]``. Then ``X_0 = B \times U``
# and ``u \in U``.
# These settings are taken from [1].
# ```math
#   \dot{x} = A \cdot x + B \cdot u
# ```

function ISSF01()
    @load ISS_path A B

    U = Hyperrectangle(low=[0.0, 0.8, 0.9], high=[0.1, 1., 1.]);  # input set
    X0 = BallInf(zeros(size(A, 1)), 0.0001)  # -0.0001 <= xi <= 0.0001 for all i
    prob_ISSF01 = @ivp(x' = A*x + B*u, x(0) ∈ X0, u ∈ U, x ∈ Universe(270))
end

# TODO:
# ISSC01 description
# The initial set ``X_0`` is defined by ``B = {x_i \in [-0001,+0.0001]``, ``\forall i}`` and
# ``U \in [0, 0.1] \times [0.8, 1.] \times [0.9, 1.]``. Then ``X_0 = B \times U``.
# These settings are taken from [1].
# ```math
#   \dot{x} = A_{ext} \cdot x
# ```

function ISSC01()
    @load ISS_path A B

    Aext = add_dimension(A, 3)
    Aext[1:270, 271:273] = B
    S = LinearContinuousSystem(Aext)

    U = Hyperrectangle(low=[0.0, 0.8, 0.9], high=[0.1, 1., 1.])  # input set
    X0 = BallInf(zeros(size(A, 1)), 0.0001)  # -0.0001 <= xi <= 0.0001 for all i
    X0 = X0 * U
    prob_ISSC01 = InitialValueProblem(S, X0)
end


# ## Results

# TODO: We use `LGG09` algorithm with sixth-order expansion in time and second order expansion
# in the spatial variables.

dirs = CustomDirections([C3, -C3]); #!jl
prob_ISSF01 = ISSF01(); #!jl
sol_ISSF01 = solve(prob_ISSF01, T=20.0, alg=LGG09(δ=6e-4, template=dirs, sparse=true, cache=false)); #!jl
#property = (-ρ(-C3, sol_ISSF01) >= -0.0007) && (ρ(C3, sol_ISSF01) <= 0.0007)

#-

using Plots, Plots.PlotMeasures, LaTeXStrings #!jl

out = [Interval(tspan(sol_ISSF01[i])) × Interval(-ρ(-C3, sol_ISSF01[i]), ρ(C3, sol_ISSF01[i])) for i in 1:10:length(sol_ISSF01)]; #!jl

fig = Plots.plot(); #!jl
Plots.plot!(fig, out, linecolor=:blue, color=:blue, alpha=0.8,
     tickfont=font(30, "Times"), guidefontsize=45,
     xlab=L"t",
     ylab=L"y_{3}",
     xtick=[0, 5, 10, 15, 20.], ytick=[-0.00075, -0.0005, -0.00025, 0, 0.00025, 0.0005, 0.00075],
     xlims=(0., 20.), ylims=(-0.00075, 0.00075),
     bottom_margin=6mm, left_margin=2mm, right_margin=4mm, top_margin=3mm,
     size=(1000, 1000)); #!jl
#nb fig

# TODO: We observe that the system converges to the equilibrium point `(1.0, 1.5)`.

dirs = CustomDirections([C3_ext, -C3_ext]); #!jl
prob_ISSC01 = ISSC01(); #!jl
sol_ISSC01 = solve(prob_ISSC01, T=20.0, alg=LGG09(δ=0.01, template=dirs, sparse=true, cache=false)); #!jl

#-

out = [Interval(tspan(sol_ISSC01[i])) × Interval(-ρ(-C3_ext, sol_ISSC01[i]), ρ(C3_ext, sol_ISSC01[i])) for i in 1:1:length(sol_ISSC01)]; #!jl

fig = Plots.plot(); #!jl
Plots.plot!(fig, out, linecolor=:blue, color=:blue, alpha=0.8, lw=1.0,
    tickfont=font(30, "Times"), guidefontsize=45,
    xlab=L"t",
    ylab=L"y_{3}",
    xtick=[0, 5, 10, 15, 20.], ytick=[-0.0002, -0.0001, 0.0, 0.0001, 0.0002],
    xlims=(0., 20.), ylims=(-0.0002, 0.0002),
    bottom_margin=6mm, left_margin=2mm, right_margin=4mm, top_margin=3mm,
    size=(1000, 1000)); #!jl
#nb fig

# TODO: We observe that the system converges to the equilibrium point `(1.0, 1.5)`.
