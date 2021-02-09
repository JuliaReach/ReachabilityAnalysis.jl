using ReachabilityAnalysis, JLD2
using ReachabilityAnalysis: add_dimension

LazySets.set_ztol(Float64, 1e-15);

examples_dir = normpath(@__DIR__, "..", "..", "..", "examples")
ISS_path = joinpath(examples_dir, "ISS", "ISS.jld2")

@load ISS_path C;
const C3 = C[3, :]; # variable y₃
const C3_ext = vcat(C3, fill(0.0, 3));

function ISSF01()
    @load ISS_path A B

    U = Hyperrectangle(low=[0.0, 0.8, 0.9], high=[0.1, 1., 1.])
    X0 = BallInf(zeros(size(A, 1)), 0.0001)
    return @ivp(x' = A*x + B*u, x(0) ∈ X0, u ∈ U, x ∈ Universe(270))
end

function ISSC01()
    @load ISS_path A B

    A_ext = add_dimension(A, 3)
    A_ext[1:270, 271:273] = B

    U = Hyperrectangle(low=[0.0, 0.8, 0.9], high=[0.1, 1., 1.])
    X0 = BallInf(zeros(size(A, 1)), 0.0001)
    X0 = X0 * U
    return @ivp(x' = A_ext*x, x(0) ∈ X0)
end

dirs = CustomDirections([C3, -C3]);
prob_ISSF01 = ISSF01();
sol_ISSF01 = solve(prob_ISSF01, T=20.0, alg=LGG09(δ=6e-4, template=dirs, sparse=true, cache=false));

πsol_ISSF01 = project(sol_ISSF01, C3);


using Plots, Plots.PlotMeasures, LaTeXStrings

fig = Plots.plot();
Plots.plot!(fig, πsol_ISSF01[1:10:end], vars=(0, 1), linecolor=:blue, color=:blue, alpha=0.8,
     xlab=L"t",
     ylab=L"y_{3}",
     xtick=[0, 5, 10, 15, 20.], ytick=[-0.00075, -0.0005, -0.00025, 0, 0.00025, 0.0005, 0.00075],
     xlims=(0., 20.), ylims=(-0.00075, 0.00075),
     bottom_margin=6mm, left_margin=2mm, right_margin=4mm, top_margin=3mm);
fig

dirs = CustomDirections([C3_ext, -C3_ext]);
prob_ISSC01 = ISSC01();
sol_ISSC01 = solve(prob_ISSC01, T=20.0, alg=LGG09(δ=0.01, template=dirs, sparse=true, cache=false));

πsol_ISSC01 = project(sol_ISSC01, C3_ext);

fig = Plots.plot();
Plots.plot!(fig, πsol_ISSC01, vars=(0, 1), linecolor=:blue, color=:blue, alpha=0.8, lw=1.0,
   #tickfont=font(30, "Times"), guidefontsize=45, !jl
    xlab=L"t",
    ylab=L"y_{3}",
    xtick=[0, 5, 10, 15, 20.], ytick=[-0.0002, -0.0001, 0.0, 0.0001, 0.0002],
    xlims=(0., 20.), ylims=(-0.0002, 0.0002),
    bottom_margin=6mm, left_margin=2mm, right_margin=4mm, top_margin=3mm);
fig

