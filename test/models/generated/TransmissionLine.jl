using ReachabilityAnalysis, LinearAlgebra, SparseArrays

function tline(;η=3, R=1.00, Rd=10.0, L=1e-10, C=1e-13 * 4.00)
    A₁₁ = zeros(η, η)
    A₁₂ = Bidiagonal(fill(-1/C, η), fill(1/C, η-1), :U)
    A₂₁ = Bidiagonal(fill(1/L, η), fill(-1/L, η-1), :L)
    A₂₂ = Diagonal(vcat(-Rd/L, fill(-R/L, η-1)))
    A  = [A₁₁ A₁₂; A₂₁ A₂₂]
    B = sparse([η+1], [1], 1/L, 2η, 1)
    return A, B
end

using Plots

A, _ = tline(η=20)
spy(A, legend=nothing, markersize=2.0, title="Sparsity pattern of A", xlab="columns", ylab="rows")

function scale!(s, α=1.0)
    s.A .*= α
    s.B .*= α
    return s
end

η = 20 # order
n = 2η # state-space dimension
A, B = tline(η=η);

Uin_ss = Interval(-0.2, 0.2)
□(ϵ) = BallInf(zeros(n), ϵ)
X0 = -inv(Matrix(A)) * B * Uin_ss ⊕ □(0.001);

Uin = Interval(0.99, 1.01)
s = @system(x'= A*x + B*u, x ∈ Universe(n), u ∈ Uin)
α = 1e-9 # scaling factor
scale!(s, α);

T = 0.7 * 1e-11 # time horizon
P = InitialValueProblem(s, X0);

sol = solve(P, T=0.7, alg=BOX(δ=1e-3));

Uout_vs_t = @. (-1.0) * project(sol, η);

fig = plot(Uout_vs_t, vars=(0, η), color=:blue, xlab="t", ylab="Uout", alpha=.5, lw=0.5)

