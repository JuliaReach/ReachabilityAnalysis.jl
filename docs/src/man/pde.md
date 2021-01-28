```@meta
DocTestSetup  = quote
    using MyPackage
end
CurrentModule = ReachabilityAnalysis
```

# Semidiscrete PDEs

## Introduction

Let's revisit linear systems but a from a different angle. Consider a system of linear differential equations of second order,

```math
    Mx''(t) + Cx'(t) + Kx(t) = f(t),\qquad t \in [0, T], \qquad (1)
```
where ``M``, ``C`` and ``K`` will be called the mass, viscosity (or damping) and stiffness matrices respectively, and ``f(t)`` is the vector of externally applied loads. A noteworthy application of such systems is the equation of equilibrium governing the linear dynamic response of a system of finite elements (FEM). Here, ``x(t)``, ``v(t) = x'(t)`` and ``a(t) = x''(t)`` are the displacement, velocity, and acceleration vectors of the finite element assemblage. Moreover, due to physical considerations it is often the case that ``M``, ``C`` and ``K`` are symmetric; ``M`` is positive-definite and ``C`` and ``K`` are positive-semidefinite. We refer to [^BAT06] for further details on the FEM context.

The initial-value problem for Eq. (1) is consists of finding a displacement ``x(t)`` satisfying Eq. (1) and given initial data,

```math
x(0) \in X_0,\qquad v(0) ∈ V_0.
```

In the rest of this section we formulate and solve problem (1) in dense time for uncertain initial conditions and input forces using linear reachability methods. We assume that ``M`` is invertible, in which case Eq. (1) can be transformed into a system of first order ODEs by multiplying with ``M^{-1}`` on the left and introducing the auxiliary vector ``[x(t),~v(t)]``.

## Free oscillations

Our first example is a damped oscillating system without forcing term,
```math
    x''(t) + 0.5~x'(t) + 4x(t) = 0, \qquad x(0) ∈ 0.7 .. 1.3, v(0) ∈ 0.2 .. 0.8.
```
The implementation and solution is straightforward. The system type [`SecondOrderContinuousSystem`](https://juliareach.github.io/MathematicalSystems.jl/latest/lib/types/#MathematicalSystems.SecondOrderLinearContinuousSystem) receives matrices
`M`, `C` and `K`; the system is assumed to be homogeneous.

```@example second_order_damped
using ReachabilityAnalysis, Plots

# x'' + 0.5x' + 4x = 0
sys = SecondOrderLinearContinuousSystem(hcat([1.0]), hcat([0.5]), hcat([4.0]))

X0 = Interval(0.7, 1.3)
V0 = Interval(0.2, 0.8)
prob = @ivp(sys, x(0) ∈ (X0, V0))
sol = solve(prob, tspan=(0.0, 20.0), alg=GLGM06(δ=0.01))

plot(sol, vars=(0, 1), lw=.1, xlab="time", lab="x(t)")
plot!(sol, vars=(0, 2), lw=0.1, xlab="time", lab="v(t)")
```
We used the zonotope-based algorithm [`GLGM06`](@ref) with step-size`δ=0.01` and
plotted the flowpipe for position $x(t)$ and velocity $v(t)$ variables.
The time horizon was chosen large enough as to show the damping effect.

As a simple illustration of working with solutions, suppose that we are interested
in the maximum and minimum values of the velocity ``v(t)`` of this system for
``t ≥ 10s``. We can filter the solution for the time interval between ``10s`` and
``20s`` with the command ``sol(10 .. 20)``. The maximum value of the velocity
is obtained by computing the support function of the flowpipe along direction
``[0.0, 1.0]``; a similar argument works for the minimum of the solution, but
with a ``-`` sign in front, as shown below.

```@example second_order_damped
vmax = ρ([0.0, 1.0], sol(10 .. 20))
vmin = -ρ([0.0, -1.0], sol(10 .. 20))

P = HPolyhedron([HalfSpace([0.0, 1.0], vmax),    # v <= vmax
                 HalfSpace([0.0, -1.0], -vmin),  # v >= vmin
                 HalfSpace([-1.0, 0.0], -10.0)]) # t >= 10

plot!(x -> 0.25, xlims=(5.0, 20.0), ylims=(-0.5, 0.8), c=:red, linestyle=:dash, lw=2.0, lab="")
plot!(x -> -0.25, xlims=(5.0, 20.0), ylims=(-0.5, 0.8), c=:red, linestyle=:dash, lw=2.0, lab="")
plot!(P, xlims=(5.0, 20.0), ylims=(-0.5, 0.8), c=:lightgreen)
```
Where we have plotted the polyhedron with max and min bounds of the flowpipe.
We have thus *proved* that ``-0.25 ≤ v(t) ≤ 0.25`` for all ``t ∈ [10, 20]`` and
for all possible values of ``x(0) ∈ X_0``, ``x'(0) ∈ V_0``.

!!! tip "Extension tip"
    To be mathematically rigorous, we would have to use set types with coefficients
    which are intervals (instead of floating-point numbers) and use the
    `IntervalArithmetic.jl` package to perform absolutely all the intermediate computations.

Let us now analyze the behavior of the solution with the damping coefficient.
If the damping coefficient is higher, then the oscillations decrease more rapidly.
Conversely, by decreasig the damping coefficient we can expect that eventually,
the property ``-0.25 ≤ v(t) ≤ 0.25`` for ``t ∈ [10, 20]`` is no longer satisfied.

```@example second_order_damped
# x'' + cx' + 4x = 0
sys_param(c) = SecondOrderLinearContinuousSystem(hcat([1.0]), hcat([c]), hcat([4.0]))
X0 = Interval(0.7, 1.3)
V0 = Interval(0.2, 0.8)
ivp(c) = @ivp(sys_param(c), x(0) ∈ (X0, V0))
sol_param(c, δ) = solve(ivp(c), tspan=(0.0, 20.0), alg=GLGM06(δ=δ));

cmin, cmax = 0.2, 0.5
cvals = range(cmin, cmax, step=0.01)
sol_tot2(δ) = sol_param.(cvals, δ);
max_v2(δ) = [ρ([0.0, 1.0], fp(10 .. 20)) for fp in sol_tot2(δ)]
min_v2(δ) = [-ρ([0.0, -1.0], fp(10 .. 20)) for fp in sol_tot2(δ)]

fig = plot(xlab="c", ylab="Velocity bounds")

plot!(fig, cvals, max_v2(0.2), lab="v max (δ=0.2)", lw=2.0, c=:green, linestyle=:dashdot)
plot!(fig, cvals, min_v2(0.2), lab="v min (δ=0.2)", lw=2.0, c=:blue, linestyle=:dashdot)

plot!(fig, cvals, max_v2(0.01), lab="v max (δ=0.01)", lw=2.0, c=:green)
plot!(fig, cvals, min_v2(0.01), lab="v min (δ=0.01)", lw=2.0, c=:blue)

plot!(fig, x -> 0.25, xlims=(cmin, cmax), ylims=(-0.5, 0.8), c=:red, linestyle=:dash, lw=2.0, lab="")
plot!(fig, x -> -0.25, xlims=(cmin, cmax), ylims=(-0.5, 0.8), c=:red, linestyle=:dash, lw=2.0, lab="")
```

## Constant loads

Consider a simple system for which the governing equations are

```math
    Mx''(t) + Kx'(t) = R,
```
where the mass matrix ``M``, stiffness matrix ``K`` and forcing term ``R`` are defined, respectively, as:
```math
M = \begin{pmatrix} 2 & 0 \\ 0 & 1 \end{pmatrix}, \qquad M = \begin{pmatrix} 6 & -2 \\ -2 & 4 \end{pmatrix},\qquad R = \begin{pmatrix} 0 \\ 10 \end{pmatrix}
```
This example is taken from Chapter 9 in [^BAT06], where the anaytic solution for
null initial conditions, ``x(0) = x'(0) = 0``, is shown to be

```math
\begin{pmatrix} x(t) \\ x'(t) \end{pmatrix} =
\begin{pmatrix} \frac{1}{\sqrt{3}} & \frac{1}{2}\sqrt{\frac{2}{3}} \\  \frac{1}{\sqrt{3}} & -\sqrt{\frac{2}{3}} \end{pmatrix}\begin{pmatrix} \frac{5}{\sqrt{3}}(1 - \cos(t\sqrt{2}) \\ 2\sqrt{\frac{2}{3}}(-1 + \cos(t\sqrt{5}) \end{pmatrix}
```
Let us first compare the analytic solution with the reachability computation
in discrete time (i.e. without bloating) and singleton initial conditions, obtained
using the `ORBIT` algorithm. This algorithm returns a sequence of singletons that
match the exact solution of the ODE at multiples of the step size. In the same
plot we show the reachability computation with step size ``δ=0.01`` using the
algorithm `GLGM06`.

```@example second_order_damped
using ReachabilityAnalysis, Plots

# Mx'' + Kx = R
M = [2 0; 0 1.]
K = [6 -2; -2 4.]
C = zeros(2, 2)
R = [0, 10.]
sys = SecondOrderAffineContinuousSystem(M, C, K, R)
X0 = zeros(2)
V0 = zeros(2)
prob = @ivp(sys, x(0) ∈ (X0, V0))

# analytic solution
A = [1/√3  (1/2)*√(2/3);
     1/√3      -√(2/3)]
x₁(t) = (5 / √3) * (1 - cos(t*√2))
x₂(t) = (2 * √(2/3)) * (-1 + cos(t*√5))
U(t) = A * [x₁(t), x₂(t)]

fig = plot(xlab="time", ylab="x1(t)")

# solution without bloating
δ = 0.1
sol_orbit = solve(prob, tspan=(0.0, 20.0), alg=ORBIT(δ=0.1))
plot!(fig, sol_orbit, vars=(0, 1), lw=2.0, markershape=:star8)

# solution with bloating
sol_lgg = solve(prob, tspan=(0.0, 20.0), alg=GLGM06(δ=0.01))
plot!(fig, sol_lgg, vars=(0, 1), lw=0.0, lab="GLGM06", c=:green)

# analytic solution
tdom = range(0.0, 20.0, step=0.001)
plot!(fig, tdom, [U(ti)[1] for ti in tdom], lab="Analytic", c=:magenta)

fig
```

The accumulated error is (numerically speaking) zero for the solution in discrete
time.


```@example second_order_damped
# check total error
x = k -> U((k-1) * δ)[1]
x′ = k -> U((k-1) * δ)[2];
```

```@example second_order_damped
sum((set(R).element[1] - x(k))^2 for (k, R) in enumerate(sol_orbit))
```

```@example second_order_damped
sum((set(R).element[2] - x′(k))^2 for (k, R) in enumerate(sol_orbit))
```

Finally, let us also note that by we can solve for a set of initial
conditions, for example, a box around zero, by changing the initial-value problem:

```@example second_order_damped
fig = plot(xlab="time", ylab="x1(t)")
X0box(ε) = BallInf(zeros(4), ε)

# solution with bloating and a box of initial conditions
ε = 1.0
prob = @ivp(sys, x(0) ∈ X0box(ε))
sol = solve(prob, tspan=(0.0, 20.0), alg=GLGM06(δ=0.01))
plot!(fig, sol, vars=(0, 1), lab="ε = $ε", lw=0.0)

ε = 0.5
prob = @ivp(sys, x(0) ∈ X0box(ε))
sol = solve(prob, tspan=(0.0, 20.0), alg=GLGM06(δ=0.01))
plot!(sol, vars=(0, 1), lab="ε = $ε", lw=0.0)

ε = 0.1
prob = @ivp(sys, x(0) ∈ X0box(ε))
sol = solve(prob, tspan=(0.0, 20.0), alg=GLGM06(δ=0.01))
plot!(sol, vars=(0, 1), lab="ε = $ε", lw=0.0)

fig
```

## Nonlinear loads

In cases where the forcing term is a nonlinear elementary function,
eg. a combination of trigonometric functions, introducing auxiliary
variables can be used to represent such function.

For instance, if ``f(t) = F \sin (3t)``, let ``u := F\sin(3t)`` and
``g := F \cos(3t)``. Then, it is a simple exercise to see that the system
```math
    x''(t) + 3.5~x'(t) + 4x(t) = F\sin(3t), \qquad x(0) ∈ 0.7 .. 1.3 × 0.7 .. 1.3
```
is formally equivalent to the following linear system:

```@example second_order_damped

F = 3.0
A = [ 0    1.0     0   0;
     -4    -3.5    1   0;
      0     0      0   F;
      0     0     -1   0]
sys = @system(x' = A * x)

B0 = BallInf(ones(2), 0.3) × Singleton([0.0]) × Interval(0.8, 2.1)
sol = solve(@ivp(sys, x(0) ∈ B0), tspan=(0.0, 15.0), alg=GLGM06(δ=0.01));

plot(sol, vars=(0, 1), c=:magenta, lw=.0, xlab="time", lab="x(t)")
plot!(sol, vars=(0, 2), c=:green, lw=.0, xlab="time", lab="v(t)")
```

## References

[^BAT06]: Bathe, Klaus-Jürgen. *Finite element procedures.* Klaus-Jurgen Bathe, 2006.
