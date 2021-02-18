using ReachabilityAnalysis, JLD2

@taylorize function lotkavolterra!(dx, x, params, t)
    local α, β, γ, δ = 1.5, 1., 3., 1.
    dx[1] = α * x[1] - β * x[1] * x[2]
    dx[2] = δ * x[1] * x[2] - γ * x[2]
    return dx
end

X0 = Hyperrectangle(low=[4.8, 1.8], high=[5.2, 2.2]);
prob = @ivp(x' = lotkavolterra!(x), dim=2, x(0) ∈ X0)
sol = solve(prob, T=8.0, alg=TMJets());
solz = overapproximate(sol, Zonotope);

using Plots

plot(solz, vars=(1, 2), alpha=0.3,lw=0., xlab="x", ylab="y", label="Flowpipe", legend=:bottomright)
plot!(X0, label="X(0)")

@taylorize function f(du, u, p, t)
    du[1] = u[3] * u[1] - u[4] * (u[1] * u[2]) - u[7] * u[1]^2
    du[2] = -u[5] * u[2] + u[6] * (u[1] * u[2])

    #encode uncertain params
    du[3] = zero(u[1]) # p[1]
    du[4] = zero(u[1]) # p[2]
    du[5] = zero(u[1]) # p[3]
    du[6] = zero(u[1]) # p[4]
    du[7] = zero(u[1]) # p[5]
end

p_int = (0.99..1.01) × (0.99..1.01) × (2.99..3.01) × (0.99..1.01) × (0.099..0.101)
U0 = Singleton([1.0, 1.0]) × convert(Hyperrectangle, p_int)
prob = @ivp(u' = f(u), dim: 7, u(0) ∈ U0);

sol = solve(prob, tspan=(0.0, 10.0));
solz = overapproximate(sol, Zonotope);

plot(solz, vars=(1, 2), lw=0.3, title="Uncertain params", lab="abs_tol = 1e-15", xlab="u1", ylab="u2")

u0 = Singleton([1.0, 1.0])
□(ϵ) = BallInf(zeros(2), ϵ)
U0 = (u0 ⊕ □(0.05)) × convert(Hyperrectangle, p_int)

prob = @ivp(u' = f(u), dim: 7, u(0) ∈ U0)

sol = solve(prob, tspan=(0.0, 10.0), TMJets(abs_tol=1e-10))
solz = overapproximate(sol, Zonotope)
plot(solz, vars=(1, 2), color=:orange, lw=0.3,
     lab="eps = 0.05", title="Uncertain u0 and uncertain params",
     xlab="u1", ylab="u2")

U0 = (u0 ⊕ □(0.01)) × convert(Hyperrectangle, p_int)
prob = @ivp(u' = f(u), dim: 7, u(0) ∈ U0)

sol = solve(prob, tspan=(0.0, 10.0), TMJets(abs_tol=1e-10))
solz = overapproximate(sol, Zonotope)
plot!(solz, vars=(1, 2), color=:blue, lw=0.3,
  lab="eps = 0.01", title="Uncertain u0 and uncertain params",
  xlab="u1", ylab="u2")

