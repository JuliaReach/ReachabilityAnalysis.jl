using ReachabilityAnalysis, Symbolics, Plots

LazySets.set_ztol(Float64, 1e-15)

function multistable_oscillator(; X0 = Interval(0.0, 0.05),
                                  V₊ = +13.5, V₋ = -13.5,
                                  R = 20.E3, C = 5.5556E-8,
                                  R1 = 20.E3, R2 = 20.E3)

    @variables x
    τ = 1/(R*C)
    α = R2/(R1+R2)

    A = -τ
    b = (τ/α) * V₊
    I₊ = HalfSpace(x <= α*V₊)
    m1 = @system(x' = Ax + b, x ∈ I₊)

    b = (τ/α) * V₋
    I₋ = HalfSpace(x >= α*V₋)
    m2 = @system(x' = Ax + b, x ∈ I₋)

    automaton = LightAutomaton(2)
    add_transition!(automaton, 1, 2, 1)
    g1 = Hyperplane(x == α*V₊)
    r1 = ConstrainedIdentityMap(1, g1)

    add_transition!(automaton, 2, 1, 2)
    g2 = Hyperplane(x == α*V₋)
    r2 = ConstrainedIdentityMap(1, g2)

    modes = [m1, m2]
    resetmaps = [r1, r2]
    H = HybridSystem(automaton, modes, resetmaps)

    # initial condition in mode 1
    X0e = [(1, X0)]
    return IVP(H, X0e)
end

prob = multistable_oscillator()

sol = solve(prob, T=100e-4, alg=INT(δ=1.E-6), fixpoint_check=false);

plot(sol, vars=(0, 1), xlab="t", ylab="v-")

tspan.(sol)

location.(sol)

plot(sol[1][end-10:end], vars=(0, 1), xlab="t", ylab="v-")
plot!(x -> 6.75, xlims=(3.1e-4, 3.3e-4), lab="Guard", lw=2.0, color=:red)

Xc = cluster(sol[1], [318, 319, 320], BoxClustering(1))

plot(sol[1][end-10:end], vars=(0, 1))
plot!(sol[2][1:10], vars=(0, 1))
plot!(x -> 6.75, xlims=(3.1e-4, 3.3e-4), lab="Guard", lw=2.0, color=:red)
plot!(Xc[1], vars=(0, 1), c=:grey)

sol = solve(prob, T=100e-4, alg=INT(δ=1.E-6), fixpoint_check=true);

plot(sol, vars=(0, 1), xlab="t", ylab="v-")

tspan(sol)

