```@meta
DocTestSetup = :(using ReachabilityAnalysis)
CurrentModule = ReachabilityAnalysis
```

# Commonly asked questions (FAQ)

### What are good introductory papers on the subject?

An elementary introduction to the principles of set-based numerical integration
can be found in [Oded Maler's](http://www-verimag.imag.fr/~maler/) article
[Computing Reachable Sets: An Introduction](https://www.semanticscholar.org/paper/Computing-Reachable-Sets-%3A-An-Introduction-Maler/299949aef669b547a36c091b768cade091d35532). For an introduction to hybrid systems reachability
we recommend the lecture notes of [Prof. Goran Frehse](https://sites.google.com/site/frehseg/),
[Formal Verification of Piecewise Affine Hybrid Systems](https://sites.google.com/site/frehseg/home#h.p_ID_132)
(DigiCosme Spring School, Paris, May 2016). Most up-to-date material related to reachability
analysis can be found in journals, conference articles or in PhD theses. We refer to
the [References](@ref all_ref) section of this manual for further links to the relevant literature.

### What happens if you consider a chaotic system?

### Can reachability analysis be used to solve large problems?

### Why did you choose Julia to write this library?

### Are there other tools that perform reachability analysis?

The wiki [Related Tools](https://github.com/JuliaReach/ReachabilityAnalysis.jl/wiki/Related-Tools)
contains an extensive list of pointers related to reachability analysis tools for
dynamical systems. A subset of such tools has participated in recent editions of
of the Friendly Competition for Applied Reachability of Continuous and Hybrid Systems,
[ARCH-COMP](https://cps-vo.org/group/ARCH).
In alphabetic order:
[CORA](https://tumcps.github.io/CORA/),
[C2E2](http://publish.illinois.edu/c2e2-tool/),
[HyDRA](https://ths.rwth-aachen.de/research/projects/hypro/),
[Hylaa](http://stanleybak.com/hylaa/),
[SpaceEx](http://spaceex.imag.fr/),
[XSpeed](http://xspeed.nitmeghalaya.in/),
[Ariadne](http://www.ariadne-cps.org/),
[DynIbex](https://perso.ensta-paris.fr/~chapoutot/dynibex/),
[Flow*](https://flowstar.org/),
[Isabelle/HOL-ODE-Numerics](https://home.in.tum.de/~hoelzl/documents/immler2012ode.pdf).

Languages and tools for hybrid systems design (as of 2006) are described in the
review article [[CPPSV06]](@ref).

### Can I use ODE solvers with interval initial conditions?

Although it is in principle possible to  ODE solvers for

### What is the wrapping effect?

Quoting a famous paper by R. E. Moore [[M65]](@ref):

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
initial set (singleton means a set with one element). For example, here we plot the
free vibration solution of a standard single degree of freedom system without physical damping,

```math
    x''(t) + 4x(t) = 0, \qquad x(0) = 1,\qquad x'(0) = 0.
```
In this initial-value problem, the initial condition is given as a point that we can
model as `X0 = Singleton([1.0, 0.0])`, where we associate the first coordinate to
position, `x(t)`, and the second coordinate to velocity, `x'(t)`.

!!! note
    Usual Julia vectors such as `X0 = [1.0, 0.0]` are also valid input, and are
    treated as a singleton. It is also valid to use tuples in second order systems,
    e.g. `prob = @ivp(X' = AX, X(0) ∈ ([1.0], [0.0]))`.

Below we plot the flowpipe for the same initial condition and different step
sizes.

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

As it is seen in the question *Can I solve a for a single initial condition?*,
even if the initial condition is a singleton, the obtained flowpipe is a sequence
of boxes in the `x-t` plane, i.e. we obtain sets with non-zero width both in time
and in space. This behavior may seem confusing at first, because the initial conditions
where determinitic. The catch is that reach-sets represents a set of states
reachable over a *time interval*, that certainly contains the exact solution for
the time-span associated to the reach-set, `tspan(R)`. The projection of the flowpipe
on the time axis thus returns a sequence of intervals, their width being the
step size of the method (in case the method has fixed step size). When we take the
Cartesian product of each time span with the projection of the flowpipe in `x(t)`,
we obtain a box.

If we consider two different step-sizes, the area of the boxes shrinks. It is known
theoretically that the flowpipe converges, in Hausdorff norm, to the exact flowpipe.
The plot below illustrates the convergence for aspects for two different step sizes,
`ΔT=0.1` and `ΔT=0.05`, evaluating the solution around the time point `3.0`.
The projection onto `x(t)` (vertical axis) shows that dividing the step size by half,
we can more accurately know the exact value of the solution, and the width of the
boxes intersecting the time point `3.0` decrease by a factor 2.5x.

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


### Does reachability solve for the vertices of the set?


### Can I compute solutions using parallel programming?

Yes. You can compute multiple flowpipes in parallel by defining an initial-value
problem with an array of initial conditions. This methods uses Julia's multithreaded parallelism,
so you have to set the number of threads to use before starting Julia.
The following example illustrates this point. For further details we refer to the section
[Distributed computations](@ref).

```@example parallel
using ReachabilityAnalysis, Plots

A = [0.0 1.0; -1.0 0.0]

B = [BallInf([0,0.] .+ k, 0.1) for k in 1:5]
prob = @ivp(x' = Ax, x(0) ∈ B)

# multi-threaded solve
sol = solve(prob, T=12.0, alg=GLGM06(δ=0.02));
plot(sol, vars=(0, 2), c=:red, alpha=.5, lw=0.2, xlab="t", ylab="y")
```
On the other hand, please note that in the example of above, you can compute with
a single integration the flowpipe corresponding to the convex hull of the elements
in the array `B`.

```@example parallel
prob = @ivp(x' = Ax, x(0) ∈ ConvexHullArray(B))
sol = solve(prob, T=12.0, alg=GLGM06(δ=0.02));
plot!(sol, vars=(0, 2), c=:lightgreen, alpha=.5, lw=0.2, xlab="t", ylab="y")
```

### How do I use the `@taylorize` macro?

The section [Some common gotchas](@ref) of the user manual details do's and dont's
for the `@taylorize` macro to speedup reachability computations using Taylor models.

### A note on interval types
