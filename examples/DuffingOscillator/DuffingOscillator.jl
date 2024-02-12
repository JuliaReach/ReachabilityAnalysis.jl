# # Duffing oscillator

#md # !!! note "Overview"
#md #     System type: Polynomial continuous system\
#md #     State dimension: 2\
#md #     Application domain: Nonlinear physics

# ## Model description

using ReachabilityAnalysis  #!jl

const ω = 1.2

@taylorize function duffing!(du, u, p, t)
    local α = -1.0
    local β = 1.0
    local δ = 0.3
    local γ = 0.37

    x, v = u

    f = γ * cos(ω * t)

    du[1] = u[2]
    du[2] = -α * x - δ * v - β * x^3 + f
    return du
end

# ## Specification

X0 = Hyperrectangle([1.0, 0.0], [0.1, 0.1])
prob = @ivp(x' = duffing!(x), x(0) ∈ X0, dim:2)

T = 2 * pi / ω;

# ## Analysis

sol = solve(prob; tspan=(0.0, 20 * T), alg=TMJets21a());

# ## Results

using Plots  #!jl
#!jl import DisplayAs  #hide

fig = plot(sol; vars=(0, 1), xlab="t", ylab="x", lw=0.2, color=:blue)

#!jl DisplayAs.Text(DisplayAs.PNG(fig))  #hide

#-

fig = plot(sol; vars=(0, 2), xlab="t", ylab="v", lw=0.2, color=:blue)

#!jl DisplayAs.Text(DisplayAs.PNG(fig))  #hide

#-

fig = plot(sol; vars=(1, 2), xlab="x", ylab="v", lw=0.5, color=:red)

#!jl DisplayAs.Text(DisplayAs.PNG(fig))  #hide
