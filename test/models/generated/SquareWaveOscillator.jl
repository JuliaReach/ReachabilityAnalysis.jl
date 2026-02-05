using ReachabilityAnalysis, Symbolics
import ReachabilityAnalysis.ReachabilityBase.Comparison as CMP

function multistable_oscillator(; X0=Interval(0.0, 0.05),
                                V₊=+13.5, V₋=-13.5,
                                R=20.E3, C=5.5556E-8,
                                R1=20.E3, R2=20.E3)
    @variables x
    τ = 1 / (R * C)
    α = R2 / (R1 + R2)
    A = -τ
    automaton = GraphAutomaton(2)

    b = (τ / α) * V₊
    I₊ = HalfSpace(x <= α * V₊)
    m1 = @system(x' = A * x + b, x ∈ I₊)

    b = (τ / α) * V₋
    I₋ = HalfSpace(x >= α * V₋)
    m2 = @system(x' = A * x + b, x ∈ I₋)

    add_transition!(automaton, 1, 2, 1)
    g1 = Hyperplane(x == α * V₊)
    r1 = ConstrainedIdentityMap(1, g1)

    add_transition!(automaton, 2, 1, 2)
    g2 = Hyperplane(x == α * V₋)
    r2 = ConstrainedIdentityMap(1, g2)

    modes = [m1, m2]
    resetmaps = [r1, r2]
    H = HybridSystem(automaton, modes, resetmaps)

    # initial condition in mode 1
    X0e = [(1, X0)]
    return IVP(H, X0e)
end;

prob = multistable_oscillator();

sol = solve(prob; T=100e-4, alg=INT(; δ=1.E-6), fixpoint_check=false);

Xc = cluster(sol[1], [318, 319, 320], BoxClustering(1));

sol = solve(prob; T=100e-4, alg=INT(; δ=1.E-6), fixpoint_check=true)
