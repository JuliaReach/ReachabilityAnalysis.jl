# # International Space Station
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/models/ISS.ipynb)
#
#md # !!! note "Overview"
#md #     System type: linear continuous system\
#md #     State dimension: 270\
#md #     Application domain: aerospace engineering
#
# ## Model description
#
# The International Space Station (ISS) is a continuous linear time-invariant
# system ``\dot{x}(t) = Ax(t) + Bu(t)`` proposed as a benchmark in ARCH
# 2016 [TLT16]. It has 270 state variables.

# The A, B, and C matrices are available in MATLAB format
# ![](slicot.org/objects/software/shared/bench-data/iss.zip)

using ReachabilityAnalysis, JLD2
using ReachabilityAnalysis: add_dimension

LazySets.set_ztol(Float64, 1e-15);
ISS_path = joinpath(@__DIR__, "ISS.jld2");

@load ISS_path C;
const C3 = C[3, :]; # variable y₃
const C3_ext = vcat(C3, fill(0.0, 3));

# ## Reachability settings
#
# Initially all the variables are in the range ``[-0.0001, 0.0001]``, ``u_1``
# is in ``[0, 0.1]``, ``u_2`` is in ``[0.8, 1]``, and ``u_3`` is in ``[0.9, 1]``.
# The time bound is 20.

# There are two versions of this benchmark:

# - ISSF01: The inputs can change arbitrarily over time: ``\forall t: u(t)\in \mathcal{U}``.

function ISSF01()
    @load ISS_path A B

    U = Hyperrectangle(low=[0.0, 0.8, 0.9], high=[0.1, 1., 1.]);  # input set
    X0 = BallInf(zeros(size(A, 1)), 0.0001)  # -0.0001 <= xi <= 0.0001 for all i
    prob_ISSF01 = @ivp(x' = A*x + B*u, x(0) ∈ X0, u ∈ U, x ∈ Universe(270))
end

# - ISSC01 (constant inputs): The inputs are uncertain only in their initial
# value, and constant over time: ``u(0)\in \mathcal{U}``, ``\dot u(t)= 0``.

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

# We use `LGG09` algorithm (based on support functions [LGG09])

# - ISSF01

dirs = CustomDirections([C3, -C3]);
prob_ISSF01 = ISSF01();
sol_ISSF01 = solve(prob_ISSF01, T=20.0, alg=LGG09(δ=6e-4, template=dirs, sparse=true, cache=false));

#-

using Plots, Plots.PlotMeasures, LaTeXStrings

out = project(sol_ISSF01[1:10:length(sol_ISSF01)], C3, vars=(0,1));
# - out = [Interval(tspan(sol_ISSF01[i])) × Interval(-ρ(-C3, sol_ISSF01[i]), ρ(C3, sol_ISSF01[i])) for i in 1:10:length(sol_ISSF01)];

fig = Plots.plot();
Plots.plot!(fig, out, linecolor=:blue, color=:blue, alpha=0.8,
#!nb     #tickfont=font(30, "Times"), guidefontsize=45, #!jl
     xlab=L"t",
     ylab=L"y_{3}",
     xtick=[0, 5, 10, 15, 20.], ytick=[-0.00075, -0.0005, -0.00025, 0, 0.00025, 0.0005, 0.00075],
     xlims=(0., 20.), ylims=(-0.00075, 0.00075),
     bottom_margin=6mm, left_margin=2mm, right_margin=4mm, top_margin=3mm);
#!nb     #size=(1000, 1000)); #!jl
fig

# - ISSC01

dirs = CustomDirections([C3_ext, -C3_ext]);
prob_ISSC01 = ISSC01();
sol_ISSC01 = solve(prob_ISSC01, T=20.0, alg=LGG09(δ=0.01, template=dirs, sparse=true, cache=false));

#-

out = project(sol_ISSC01, C3_ext, vars=(0,1));
# - out = [Interval(tspan(sol_ISSC01[i])) × Interval(-ρ(-C3_ext, sol_ISSC01[i]), ρ(C3_ext, sol_ISSC01[i])) for i in 1:1:length(sol_ISSC01)];

fig = Plots.plot();
Plots.plot!(fig, out, linecolor=:blue, color=:blue, alpha=0.8, lw=1.0,
#!nb    #tickfont=font(30, "Times"), guidefontsize=45, !jl
    xlab=L"t",
    ylab=L"y_{3}",
    xtick=[0, 5, 10, 15, 20.], ytick=[-0.0002, -0.0001, 0.0, 0.0001, 0.0002],
    xlims=(0., 20.), ylims=(-0.0002, 0.0002),
    bottom_margin=6mm, left_margin=2mm, right_margin=4mm, top_margin=3mm);
#!nb    #, size=(1000, 1000)); #!jl
fig

# ## References

# - [TLT16] Tran, Hoang-Dung, Luan Viet Nguyen, and Taylor T. Johnson. "Large-scale linear systems from order-reduction (benchmark proposal)." 3rd Applied Verification for Continuous and Hybrid Systems Workshop (ARCH), Vienna, Austria. 2016.
# - [LGG09] Le Guernic, Colas, and Antoine Girard. "Reachability analysis of linear systems using support functions." Nonlinear Analysis: Hybrid Systems 4.2 (2010): 250-262.
