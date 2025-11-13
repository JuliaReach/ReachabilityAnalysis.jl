using ReachabilityAnalysis

function opamp_nondet(; X0=Singleton([0.0]),
                      R₁=2.0, R₂=6.0, C=1.e-3,
                      Ein=Interval(1.9, 2.1))
    α = hcat(-1 / (R₂ * C))
    β = hcat(-1 / (R₁ * C))

    # continuous evolution
    sys = @system(eout' = α * eout + β * ein, ein ∈ Ein, eout ∈ Universe(1))

    # initial-value problem
    return @ivp(sys, x(0) ∈ X0)
end;

sol_nondet = solve(opamp_nondet(); T=0.1, alg=INT(; δ=1e-4));

sol_nondet = solve(opamp_nondet(; X0=Interval(-1.0, 1.0)); T=0.1, alg=INT(; δ=1e-4))

sol_nondet = solve(opamp_nondet(; X0=Interval(-0.5, 0.5)); T=0.1, alg=INT(; δ=1e-4))

sol_nondet = solve(opamp_nondet(); T=0.1, alg=INT(; δ=1e-4))

using Symbolics

function opamp_with_saturation(; X0=Singleton(zeros(2)),
                               R₁=2.0, R₂=6.0, C=1.e-3,
                               γ=100.0, δ=100.0, Es=2.0)
    var = @variables eₒ ein
    α = hcat(-1 / (R₂ * C))
    β = hcat(-1 / (R₁ * C))

    # transition graph
    automaton = GraphAutomaton(2)
    add_transition!(automaton, 1, 2, 1)

    # modes
    A = [α β; 0 γ]
    b = [0; δ]
    mode1 = @system(x' = A * x + b, x ∈ HalfSpace(ein <= Es, var))
    A = [α β; 0 0]
    mode2 = @system(x' = A * x, x ∈ Universe(2))
    modes = [mode1, mode2]

    # transition mode1 -> mode2 (saturation)
    t1 = @map(x -> x, dim:2, x ∈ Hyperplane(ein == Es, var))
    resetmaps = [t1]

    # initial condition: at the origin in mode 1
    initial_condition = [(1, X0)]

    H = HybridSystem(automaton, modes, resetmaps, [AutonomousSwitching()])
    return IVP(H, initial_condition)
end;

X0 = Hyperrectangle([0.0, 1.5], [0.0, 0.1]);

prob_const = opamp_with_saturation(; X0=X0, γ=0.0, δ=0.0, Es=3.0);

sol_const = solve(prob_const; T=0.1, alg=BOX(; δ=1e-4));

X0 = Singleton(zeros(2));

prob_lin = opamp_with_saturation(; X0=X0, γ=0.0, δ=100.0, Es=1.0)
prob_exp = opamp_with_saturation(; X0=X0, γ=-100.0, δ=100.0, Es=3.0);

sol_lin = solve(prob_lin; T=0.1, alg=BOX(; δ=1e-4))

sol_exp = solve(prob_exp; T=0.1, alg=BOX(; δ=1e-4));
