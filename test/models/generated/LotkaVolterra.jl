@taylorize function lotkavolterra!(du, u, p, t)
    local α, β, γ, δ = 1.5, 1.0, 3.0, 1.0

    x, y = u
    xy = x * y
    du[1] = α * x - β * xy
    du[2] = δ * xy - γ * y
    return du
end

X0 = Hyperrectangle(; low=[4.8, 1.8], high=[5.2, 2.2])
prob = @ivp(x' = lotkavolterra!(x), dim:2, x(0) ∈ X0);

sol = solve(prob; T=8.0, alg=TMJets())
solz = overapproximate(sol, Zonotope);

fig = plot(solz; vars=(1, 2), alpha=0.3, lw=0.0, xlab="x", ylab="y",
           lab="Flowpipe", legend=:bottomright)
plot!(fig, X0; label="X(0)")

@taylorize function lotkavolterra_parametric!(du, u, p, t)
    x, y, αp, βp, γp, δp, ϵp = u
    xy = x * y
    du[1] = αp * x - βp * xy - ϵp * x^2
    du[2] = δp * xy - γp * y

    # encode uncertain parameters
    du[3] = zero(αp)
    du[4] = zero(βp)
    du[5] = zero(γp)
    du[6] = zero(δp)
    du[7] = zero(ϵp)
    return du
end

p_int = (0.99 .. 1.01) × (0.99 .. 1.01) × (2.99 .. 3.01) × (0.99 .. 1.01) × (0.099 .. 0.101)
U0 = cartesian_product(Singleton([1.0, 1.0]), convert(Hyperrectangle, p_int))
prob = @ivp(u' = lotkavolterra_parametric!(u), dim:7, u(0) ∈ U0);

sol = solve(prob; tspan=(0.0, 10.0))
solz = overapproximate(sol, Zonotope);

fig = plot(solz; vars=(1, 2), lw=0.3, title="Uncertain parameters",
           lab="abstol = 1e-15", xlab="x", ylab="y")

□(ϵ) = BallInf([1.0, 1.0], ϵ)

U0 = cartesian_product(□(0.05), convert(Hyperrectangle, p_int))
prob = @ivp(u' = lotkavolterra_parametric!(u), dim:7, u(0) ∈ U0);

sol = solve(prob; T=10.0, alg=TMJets(; abstol=1e-10))
solz = overapproximate(sol, Zonotope);

fig = plot(solz; vars=(1, 2), color=:orange, lw=0.3, lab="ϵ = 0.05",
           title="Uncertain u0 and uncertain parameters", xlab="x", ylab="y")

U0 = cartesian_product(□(0.01), convert(Hyperrectangle, p_int))
prob = @ivp(u' = lotkavolterra_parametric!(u), dim:7, u(0) ∈ U0);

sol = solve(prob; T=10.0, alg=TMJets(; abstol=1e-10))
solz = overapproximate(sol, Zonotope);

plot!(solz; vars=(1, 2), color=:blue, lw=0.3, lab="ϵ = 0.01")
