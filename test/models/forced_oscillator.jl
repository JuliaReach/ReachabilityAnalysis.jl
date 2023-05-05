# See [Chapter 9, in Finite Element Procedures, K-J Bathe].
function forced_oscillator(; X0=(zeros(2), zeros(2)), X=Universe(2))
    # Mx'' + Kx = R
    M = [2 0; 0 1.0]
    K = [6 -2; -2 4.0]
    C = zeros(2, 2)
    R = [0, 10.0]
    example_9_1_Bathe = SecondOrderAffineContinuousSystem(M, C, K, R)
    prob = InitialValueProblem(example_9_1_Bathe, X0)
    tspan = (0.0, 20.0)
    return prob, tspan
end

# analytic solution (see Example 9.7 in Bathe)
function forced_oscillator_solution()
    A = [1/√3 (1 / 2)*√(2 / 3);
         1/√3 -√(2 / 3)]
    x₁(t) = (5 / √3) * (1 - cos(t * √2))
    x₂(t) = (2 * √(2 / 3)) * (-1 + cos(t * √5))
    return U(t) = A * [x₁(t), x₂(t)]
end
