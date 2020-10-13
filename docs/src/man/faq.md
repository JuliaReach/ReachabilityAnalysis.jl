# Commonly asked questions (FAQ)

### What are good introductory papers on the subject?

### What happens if you consider a chaotic system?

### Can reachability analysis be used to solve large problems?

### Are there other tools that perform reachability analysis?

Yes. The wiki [Related Tools](https://github.com/JuliaReach/ReachabilityAnalysis.jl/wiki/Related-Tools)
contains an extensive list of tools that perform reachability analysis.

### Can I use ODE solvers with interval initial conditions?

### Why did you choose Julia to write this library?

### What is the wrapping effect?

Quoting a famous paper by Moore [[M65]](@ref):

> Under the flow itself a box is carried after certain time into a set of points
> which will in general not remain a box excepted for a few simple flows.

The wrapping effect is associated to error propagation that arises from the inability
of a set representation to accurately abstract some properties of system
under study. Wrapping-free methods only exist for linear systems. Nonlinear reachability
methods control, but do not totally elimitate, the wrapping effect in different ways.

As a simple illustration of the wrapping effect, consider the image of a box
through a rotation in discrete time-steps. Box enclosures (full line)
introduce a strong wrapping effect. On the other hand, had we represented the
sequence using zonotopes (dashed lines), the result would be exact, i.e. without
wrapping, since the image of a hyperrectangular set under an affine map *is* a zonotope.

```@example
using LazySets, Plots

R(θ) = [cos(θ) -sin(θ); sin(θ) cos(θ)]

B = BallInf(ones(2), 0.1)
B′ = R(π/4) * B
B′′ = R(π/2) * B

□B = box_approximation(B)
□B′ = box_approximation(B′)
□B′′ = box_approximation(R(π/4) * □B′)

plot(B, ratio=1, lw=2.0, style=:dash)
plot!(B′, lw=2.0, style=:dash)
plot!(B′′, lw=2.0, style=:dash)

plot!(□B, lw=2.0, style=:solid)
plot!(□B′, lw=2.0, style=:solid)
plot!(□B′′, lw=2.0, style=:solid)
```

### Can I solve a for a single initial condition?

To solve for a single initial condition, i.e. a "point", use `Singleton` as the
initial set; singleton means a set with one element. For example, here we plot the
free vibration solution of a standard single degree of freedom system without physical damping,

```math
    x''(t) + 4x(t) = 0, \qquad x(0) = 1,\qquad x'(0) = 0.
```
In this initial-value problem, the initial condition is given as a point that we can
model as `Singleton([1.0, 0.0])`.

!!! note
    Vectors such as `X0 = [1.0, 0.0]` are also a valid input, and are treated as
    a singleton.

```@example cosine
using ReachabilityAnalysis, Plots

# x' = v
# v' = -4x
A = [0 1; -4. 0]
X0 = Singleton([1.0, 0.0])
prob = @ivp(X' = AX, X(0) ∈ X0)

f(ΔT) = solve(prob, tspan=(0.0, 5.0), alg=GLGM06(δ=ΔT))

plot(f(0.3), vars=(0, 1), lab="ΔT=0.3", color=:yellow)
plot!(f(0.1), vars=(0, 1), lab="ΔT=0.1", color=:lightblue)
plot!(f(0.05), vars=(0, 1), xlab="time", ylab="x(t)", lab="ΔT=0.05", color=:green)

dom = 0:0.01:5.0
plot!(dom, cos.(2.0 * dom), lab="Analytic", color=:magenta)
```

### Why do I see boxes for single initial conditions?

Even for singleton initial conditions, a reachability computation returns a
sequence of *sets*.

```@example cosine
plot(f(0.1)(3.0), vars=(0, 1), xlab="time", ylab="x(t)", lab="ΔT=0.1", color=:lightblue)

I(Δt, t) = -ρ([-1.0, 0.0], f(Δt)(t)) .. ρ([1.0, 0.0], f(Δt)(t)) |> Interval

I01 = I(0.1, 3.0)
plot!(y -> max(I01), xlims=(2.9, 3.1), lw=3.0, style=:dash, color=:lightblue, lab="Δx = $(I01.dat)")
plot!(y -> min(I01), xlims=(2.9, 3.1), lw=3.0, style=:dash, color=:lightblue, lab="")
plot!(f(0.05)(3.0), vars=(0, 1), xlab="time", ylab="x(t)", lab="ΔT=0.05", color=:green)

I005 = I(0.05, 3.0)
plot!(y -> max(I005), xlims=(2.9, 3.1), lw=3.0, style=:dash, color=:green, lab="Δx = $(I005.dat)")
plot!(y -> min(I005), xlims=(2.9, 3.1), lw=3.0, style=:dash, color=:green, lab="")

dom = 2.9:0.01:3.1
plot!(dom, cos.(2.0 * dom), lab="Analytic", color=:magenta, legend=:bottomright)
```

### Why do some trajectories escape the flowpipe?
