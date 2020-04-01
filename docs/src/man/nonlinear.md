# Nonlinear ordinary differential equations

In this section we illustrate the flowpipe computation for a nonlinear system.

## Model description

Our running example is the [Lotka-Volterra model](https://en.wikipedia.org/wiki/Lotka%E2%80%93Volterra_equations). The 2-dimensional Lotka-Volterra system depicts the populations change of a class of predators and a class of
preys. The growth rate of preys’ population $x$ over time is given by

```math
\dot{x} = x\cdot (\alpha - \beta \cdot y)
```
wherein  $\alpha, \beta$ are constant parameters and $y$ is the population of predators.

It gives that the number of preys grows exponentially without predation.

The population growth of predators is governed by the differential equation

```math
\dot{y} = -y\cdot (\gamma - \delta\cdot x)
```
wherein  $\gamma, \delta$ are constant parameters.

We set those parameters as  $\alpha = 1.5 ,  \beta = 1 ,  \gamma = 3$  and  $\delta = 1$.

```@example lotka_volterra
using ReachabilityAnalysis

@taylorize function lotka_volterra!(du, u, p, t)
    local α, β, γ, δ = 1.5, 1.0, 3.0, 1.0
    du[1] = u[1] * (α - β*u[2])
    du[2] = -u[2] * (γ - δ*u[1])
    return du
end
```

## Reachability settings

The reachability settings are taken from [this resource](https://ths.rwth-aachen.de/research/projects/hypro/lotka-volterra/).

We consider the initial set  $x\in [4.8,5.2], y \in [1.8,2.2]$.

```@example lotka_volterra
X₀ = Hyperrectangle(low=[4.8, 1.8], high=[5.2, 2.2])

prob = @ivp(x' = lotka_volterra!(x), dim: 2, x(0) ∈ X₀)
```

## Results

We compute the flowpipe using the TMJets algorithm for the time horizon $[0,5]$:

```@example lotka_volterra
sol = solve(prob, T=5.0)

setrep(sol)
```

We can change to the zonotopic overapproximation of the flowpipe using
the `overapproximate` function:

```@example lotka_volterra
sol = overapproximate(sol, Zonotope)

setrep(sol)
```

Finally we plot the solution in phase-space:

```@example lotka_volterra
using Plots

plot(sol, vars=(1, 2), xlab="x", ylab="y", lw=0.2, color=:lightblue, lab="Flowpipe")
plot!(X₀, color=:orange, lab="Xo")
```
