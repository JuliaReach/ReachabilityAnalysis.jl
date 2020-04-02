# =======================================
# Van der Pol oscillator
# #include("models/lotka_volterra.jl")
# https://en.wikipedia.org/wiki/Van_der_Pol_oscillator
# Number of state variables: 2
# =======================================

@taylorize function vanderPol!(dx, x, params, t)
    local μ = 1.0
    dx[1] = x[2]
    dx[2] = (μ * x[2]) * (1 - x[1]^2) - x[1]
    return dx
end

function vanderpol()
    X0 = Hyperrectangle(low=[1.25, 2.35], high=[1.55, 2.45])
    prob = @ivp(x' = vanderPol!(x), dim: 2, x(0) ∈ X0)
    tspan = (0.0, 5.0)

    return prob, tspan
end

#=
OLD

# check mode
property = (t, x) -> x[2] <= 2.75
sol = solve(P, T=7.0, property=property, alg)

# test set representation option
sol = solve(P, T=7.0, setrep=:Hyperrectangle, alg)
@test set(sol[1]) isa Hyperrectangle

sol = solve(P, T=7.0, setrep=:IntervalBox, alg)
@test set(sol[1]) isa IntervalBox

sol = solve(P, T=7.0, setrep=:Zonotope, alg)
@test set(sol[1]) isa Zonotope

sol = solve(P, T=7.0, setrep=:TaylorModel, alg)
@test set(sol[1]) isa TaylorModel
=#
