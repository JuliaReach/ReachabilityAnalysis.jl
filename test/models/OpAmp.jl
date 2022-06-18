function opamp_nondet(; X0 = Interval(0.0, 0.0),
                        R₁ = 2., R₂ = 6., C = 1.e-3,
                        Ein = Interval(1.9, 2.1))

    α = hcat(-1/(R₂ * C))
    β = hcat(-1/(R₁ * C))

    # continuous evolution
    invariant = Universe(1) # no state constraint
    sys = @system(eout' = α*eout + β*ein, ein ∈ Ein, eout ∈ invariant)

    # initial-value problem
    return @ivp(sys, x(0) ∈ X0)
end

sol_nondet = solve(opamp_nondet(), T=0.1, alg=INT(δ=1e-4));


fig = plot()
plot!(fig, sol_nondet, vars=(0, 1), xlab=L"t", ylab=L"e_{out}", title="Solution for non-deterministic input", lw=0.2)

fig = plot(xlab=L"t", ylab=L"e_{out}(t)", title="Solution for non-deterministic input")

Δt = 0 .. 0.04

sol = solve(opamp_nondet(X0=Interval(-1.0, 1.0)), T=0.1, alg=INT(δ=1e-4));
plot!(fig, sol(Δt), vars=(0, 1), lab="X0 = -1 .. 1", lw=0.2)

sol = solve(opamp_nondet(X0=Interval(-0.5, 0.5)), T=0.1, alg=INT(δ=1e-4));
plot!(fig, sol(Δt), vars=(0, 1), lab="X0 = -0.5 .. 0.5", lw=0.2)

sol = solve(opamp_nondet(), T=0.1, alg=INT(δ=1e-4));
plot!(fig, sol(Δt), vars=(0, 1), lab="X0 = 0", lw=0.2)

fig

using Symbolics

function opamp_with_saturation(; X0 = BallInf(zeros(2), 0.0),
                                 R₁ = 2., R₂ = 6., C = 1.e-3,
                                 γ = 100., δ = 100., Es = 2.)
    var = @variables eₒ ein
    α = hcat(-1/(R₂ * C))
    β = hcat(-1/(R₁ * C))

    # transition graph
    automaton = GraphAutomaton(2)
    add_transition!(automaton, 1, 2, 1)

    # modes
    A = [α β; 0. γ]
    b = [0.; δ]
    invariant = HalfSpace(ein <= Es, var)
    mode1 = @system(x' = Ax + b, x ∈ invariant)
    A = [α β; 0. 0.];
    mode2 = @system(x' = Ax, x ∈ Universe(2))
    modes = [mode1, mode2]

    # transition mode1 -> mode2 (saturation)
    guard = Hyperplane(ein == Es, var)
    t1 = @map(x -> x, dim:2, x ∈ guard)
    resetmaps = [t1]

    H = HybridSystem(automaton, modes, resetmaps, [AutonomousSwitching()])

    # initial condition is at the origin in mode 1
    initial_condition = [(1, X0)]
    return IVP(H, initial_condition)
end

X0 = Hyperrectangle([0., 1.5], [0.0, 0.1]);

prob_const = opamp_with_saturation(X0=X0, γ = 0., δ = 0., Es = 3.)
sol_const = solve(prob_const, T=0.1, alg=BOX(δ=1e-4));

plot(sol_const, vars=(0, 2), xlab=L"t", ylab=L"e_{in}(t)", title="Constant input signal", lw=0.2)
plot!(x->x, x->1.5, 0., 0.1, line=2, color="red", linestyle=:dash, legend=nothing)

plot(sol_const, vars=(0, 1), xlab=L"t", ylab=L"e_{out}(t)", title="Solution for constant input", lw=0.2)

X0 = Singleton(zeros(2));

# linearly increasing input signal
prob_lin = opamp_with_saturation(X0=X0, γ=0., δ=100., Es=1.);
sol_lin = solve(prob_lin, T=0.1, alg=BOX(δ=1e-4));

# exponentially increasing input signal
prob_exp = opamp_with_saturation(X0=X0, γ=-100., δ=100., Es=3.);
sol_exp = solve(prob_exp, T=0.1, alg=BOX(δ=1e-4));

plot(sol_lin, vars=(0, 2), xlab=L"t", ylab=L"e_{in}(t)", title="Linear input w/saturation", lw=0.2)

plot(sol_exp, vars=(0, 2), xlab=L"t", ylab=L"e_{in}(t)", title="Exponential input", lw=0.2)

plot(sol_lin, vars=(0, 1), xlab=L"t", ylab=L"e_{out}(t)", title="Solution for linear input", lw=0.2)

plot(sol_exp, vars=(0, 1), xlab=L"t", ylab=L"e_{out}(t)", title="Solution for exponential input", lw=0.2)

