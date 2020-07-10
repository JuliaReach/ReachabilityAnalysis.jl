#TODO:

using ReachabilityAnalysis, JLD2
using ReachabilityAnalysis: add_dimension


const C3 = C[3, :]; #!jl # variable y₃

function ISSF01()
    @load ISS_path A B

    U = Hyperrectangle(low=[0.0, 0.8, 0.9], high=[0.1, 1., 1.]);  # input set
    X0 = BallInf(zeros(size(A, 1)), 0.0001)  # -0.0001 <= xi <= 0.0001 for all i
    prob_ISSF01 = @ivp(x' = A*x + B*u, x(0) ∈ X0, u ∈ U, x ∈ Universe(270))
end

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

#property = (-ρ(-C3, sol_ISSF01) >= -0.0007) && (ρ(C3, sol_ISSF01) <= 0.0007)

Plots.plot!(fig, out, linecolor=:blue, color=:blue, alpha=0.8,
     tickfont=font(30, "Times"), guidefontsize=45,
     xlab=L"t",
     ylab=L"y_{3}",
     xtick=[0, 5, 10, 15, 20.], ytick=[-0.00075, -0.0005, -0.00025, 0, 0.00025, 0.0005, 0.00075],
     xlims=(0., 20.), ylims=(-0.00075, 0.00075),
     bottom_margin=6mm, left_margin=2mm, right_margin=4mm, top_margin=3mm,

Plots.plot!(fig, out, linecolor=:blue, color=:blue, alpha=0.8, lw=1.0,
    tickfont=font(30, "Times"), guidefontsize=45,
    xlab=L"t",
    ylab=L"y_{3}",
    xtick=[0, 5, 10, 15, 20.], ytick=[-0.0002, -0.0001, 0.0, 0.0001, 0.0002],
    xlims=(0., 20.), ylims=(-0.0002, 0.0002),
    bottom_margin=6mm, left_margin=2mm, right_margin=4mm, top_margin=3mm,

