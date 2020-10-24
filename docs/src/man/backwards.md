# Backwards reachability

```@example
using ReachabilityAnalysis, Plots

@taylorize function rotating!(du, u, p, t)
    du[1] = u[2]
    du[2] = -u[1]
end

# A = [0. 1.0; -1 0]
# prob = @ivp(x' = Ax, x(0) ∈ BallInf(ones(2), 0.5))
prob = @ivp(x' = rotating!(x), dim: 2, x(0) ∈ BallInf(ones(2), 0.5))

sol_forward = solve(prob, tspan=(0.0, 5.0));
sol_backward = solve(prob, tspan=(5.0, 0.0));

solz_forward = overapproximate(sol_forward, Zonotope)
solz_backward = overapproximate(sol_backward, Zonotope);
```
