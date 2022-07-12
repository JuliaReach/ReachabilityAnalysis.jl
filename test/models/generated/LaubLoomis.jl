using ReachabilityAnalysis, Plots

@taylorize function laubloomis!(dx, x, params, t)
    dx[1] = 1.4*x[3] - 0.9*x[1]
    dx[2] = 2.5*x[5] - 1.5*x[2]
    dx[3] = 0.6*x[7] - 0.8*(x[2]*x[3])
    dx[4] = 2 - 1.3*(x[3]*x[4])
    dx[5] = 0.7*x[1] - (x[4]*x[5])
    dx[6] = 0.3*x[1] - 3.1*x[6]
    dx[7] = 1.8*x[6] - 1.6*(x[2]*x[7])
    return dx
end

function laubloomis(; W=0.01)
    # initial states
    X0c = [1.2, 1.05, 1.5, 2.4, 1.0, 0.1, 0.45]
    X0 = Hyperrectangle(X0c, fill(W, 7))

    # initial-value problem
    prob = @ivp(x' = laubloomis!(x), dim: 7, x(0) ∈ X0)

    return prob
end

prob = laubloomis(W=0.01)
alg = TMJets(abstol=1e-11, orderT=7, orderQ=1, adaptive=true);

sol_1 = solve(prob, T=20.0, alg=alg);

sol_1z = overapproximate(sol_1, Zonotope);

# canonical direction along x₄
const e4 = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0];

ρ(e4, sol_1z)

ρ(e4, sol_1z) < 4.5

@show ρ(e4, sol_1z[end]) + ρ(-e4, sol_1z[end])

prob = laubloomis(W=0.05)
alg = TMJets(abstol=1e-12, orderT=7, orderQ=1, adaptive=false);

sol_2 = solve(prob, T=20.0, alg=alg)
sol_2z = overapproximate(sol_2, Zonotope);

ρ(e4, sol_2z)

ρ(e4, sol_2z) < 4.5

# width of final box
@show ρ(e4, sol_2z[end]) + ρ(-e4, sol_2z[end])

prob = laubloomis(W=0.1)
alg = TMJets(abstol=1e-12, orderT=7, orderQ=1, adaptive=false);

sol_3 = solve(prob, T=20.0, alg=alg)
sol_3z = overapproximate(sol_3, Zonotope);

ρ(e4, sol_3z)

ρ(e4, sol_3z) < 5.0

# width of final box
@show ρ(e4, sol_3z[end]) + ρ(-e4, sol_3z[end])


fig = plot()
plot!(fig, sol_3z, vars=(0, 4), linecolor="green", color=:green, alpha=0.8)
plot!(fig, sol_2z, vars=(0, 4), linecolor="blue",  color=:blue, alpha=0.8)
plot!(fig, sol_1z, vars=(0, 4), linecolor="yellow", color=:yellow, alpha=0.8,
      tickfont=font(10, "Times"), guidefontsize=15,
      xlab=L"t",
      ylab=L"x_4",
      xtick=[0., 5., 10., 15., 20.], ytick=[1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5.],
      xlims=(0., 20.), ylims=(1.5, 5.02),
      bottom_margin=6mm, left_margin=2mm, right_margin=4mm, top_margin=3mm,
      size=(600, 600))

plot!(fig, x->x, x->4.5, 0., 20., line=2, color="red", linestyle=:dash, legend=nothing)
plot!(fig, x->x, x->5., 0., 20., line=2, color="red", linestyle=:dash, legend=nothing)

