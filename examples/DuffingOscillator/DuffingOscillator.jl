# # Duffing Oscillator
#
#md # !!! note "Overview"
#md #     System type: nonlinear\
#md #     State dimension: 2\
#md #     Application domain: Nonlinear physics
#
# ## Implementation

using ReachabilityAnalysis, Plots

@taylorize function duffing!(du, u, p, t)
    local α = -1.0
    local β = 1.0
    local δ = 0.3
    local γ = 0.37

    x, v = u
    f = γ * cos(ω * t)

    du[1] = u[2]
    return du[2] = -α * x - δ * v - β * x^3 + f
end

ω = 1.2
T = 2 * pi / ω
X0 = Singleton([1.0, 0.0]) ⊕ BallInf(zeros(2), 0.1)
prob = @ivp(x' = duffing!(x), x(0) ∈ X0, dim = 2);

# ## Results

sol = solve(prob; tspan=(0.0, 20 * T), alg=TMJets21a());

#-

fig = plot(sol; vars=(0, 1), xlab="t", ylab="x", lw=0.2, color=:blue)

#!jl import DisplayAs  #hide
#!jl DisplayAs.Text(DisplayAs.PNG(fig))  #hide

#-

fig = plot(sol; vars=(0, 2), xlab="t", ylab="v", lw=0.2, color=:blue)

#!jl DisplayAs.Text(DisplayAs.PNG(fig))  #hide

#-

fig = plot(sol; vars=(1, 2), xlab="x", ylab="v", lw=0.5, color=:red)

#!jl DisplayAs.Text(DisplayAs.PNG(fig))  #hide
