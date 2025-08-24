using ReachabilityAnalysis, LazySets, Plots, IntervalMatrices

const N = Float64

# System matrix with uncertainty
AS = MatrixZonotope(N[-1 -5; 1 -1], [N[0 0; 0.1 0]], [6])

# Initial set as a box, converted to SparsePolynomialZonotope
X0 = Zonotope(N[1.0, 1.0], [0.2 0; 0 0.2])
X0 = convert(SparsePolynomialZonotope, X0)

# Time step
δ = N(2π) / 200
T = 2 * N(π)

# Parameters to vary
tols = [1e-2, 1e-12]
recursives = [true, false]

# Prepare 2x2 layout
plot_layout = @layout [a b; c d]
p = plot(layout=plot_layout, size=(900,800))

# Iterate over parameter combinations
for (i_tol, tol) in enumerate(tols)
    for (i_rec, recursive) in enumerate(recursives)
        alg = HLBS25(
            δ = δ,
            approx_model = CorrectionHullMatrixZonotope(),
            max_order = 4,
            taylor_order = 10,
            reduction_method = LazySets.GIR05(),
            recursive = recursive,
            tol = tol,
            norm = Inf
        )
        prob = @ivp(x' = A * x, x(0) ∈ X0, A ∈ AS)
        sol = solve(prob, alg; T=T)
        
        # Plot every 5th step for clarity
        subplot = (i_tol-1)*2 + i_rec
        plot!(p[subplot], X0, vars=(1,2), color=:red, alpha=0.3, lab="X₀")
        for (j, rset) in enumerate(sol)
            if j % 5 == 0
                plot!(p[subplot], set(rset), vars=(1,2),
                      lw=0.5, color=:blue, alpha=0.8, lab="")
            end
        end
        title!(p[subplot], "tol=$tol, recursive=$recursive")
        xlabel!(p[subplot], "x₁")
        ylabel!(p[subplot], "x₂")
    end
end

display(p)
