using ReachabilityAnalysis, LazySets, Plots

const N = Float64

AS = MatrixZonotope(N[0 -1; 1 0], [zeros(N, 2, 2)], [5])

X0 = SparsePolynomialZonotope(N[1.0, 1.0],         # center
                              N(0.1) * N[2.0 0.0 1.0; 1.0 2.0 1.0],  # generators
                              N(0.1) * reshape(N[1.0, 0.5], 2, 1), # dependent generators
                              [1 0 1; 0 1 3])      # exponents

prob = @ivp(x' = A * x, x(0) ∈ X0, A ∈ AS)
δ = N(2π) / 200

# Instantiate the HLBS25 algorithm.
# - δ: Step size.
# - approx_model: Discretization model. CorrectionHullMatrixZonotope is suitable for matrix zonotopes.
# - max_order: Maximum order for the Taylor series expansion of the matrix exponential.
# - taylor_order: Taylor series order for each step.
# - reduction_method: Method for reducing the order of matrix zonotopes.
# - recursive: Specifies whether to use a recursive method for the Taylor expansion.
alg = HLBS25(δ=δ,
             approx_model=CorrectionHullMatrixZonotope(),
             max_order=5,
             taylor_order=5,
             reduction_method=LazySets.GIR05(),
             recursive=false)
T = N(π)
sol = solve(prob, alg; T=T)

plot(sol, vars=(1, 2), lw=0.5, color=:blue, alpha=0.8,
     xlab="x₁", ylab="x₂",
     title="Flowpipe of the parametric linear system using HLBS25",
     lab="Flowpipe", legend=:bottomright)

# Plot the initial set X0 for context.
plot!(X0, vars=(1, 2), color=:red, alpha=0.5, lab="X₀")
