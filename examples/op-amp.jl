# # International Space Station
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/models/op-amp.ipynb)
#
#md # !!! note "Overview"
#md #     System type: linear continuous system\
#md #     State dimension: 2\
#md #     Application domain: electrical engineering
#
# ## Model description
#

# Operational amplifiers ("op-amp") are electronic devices to amplify an input
# voltage signal. They are widely used to construct filters that remove a
# desired range of frequencies from the input signal [1].

# The schematic diagram of an op-amp is shown below. It has two terminals:
# - input side (left)
# - output terminal (right)

# ![](https://github.com/JuliaReach/ReachabilityAnalysis.jl/blob/master/examples/fig/opamp.png?raw=true)

# The output voltage is ``e_o = K(e_B - e_A)``, where ``K`` is the voltage gain
# of the op-amp. ``K`` is usually vey large (the order of ``10^5 V / V``).


# In the *ideal op-amp* scenario, the following assumptions hold:
# - The input terminals of the op-amp draw negligible current.
# - The voltage difference at the input terminals ``e_B - e_A`` is zero.
# - The gain ``K`` is infinite.

# ## Example (inverting amplifier)

# ![](https://github.com/JuliaReach/ReachabilityAnalysis.jl/blob/master/examples/fig/opamp_circuit.png?raw=true)

# We can show that: ``e_o = -\dfrac{KR_2}{R_1 + R_2 + KR_1}e_{in}(t)``

# In the ideal op-amp, $K\to \infty$ and ``e_O = -\dfrac{R_2}{R_1}e_{in}(t)``

# ## Case of study (inverter with a capacitor)

# ![](https://github.com/JuliaReach/ReachabilityAnalysis.jl/blob/master/examples/fig/opamp_c.png?raw=true)

# We can show that the output satisfied the following ODE

# ``\dfrac{de_o(t)}{dt} = -\dfrac{1}{R_2 C}e_o(t) - \dfrac{1}{R_1C}e_{in}(t)``

# - Equation for the input voltage:

# ``\dfrac{d \ ein}{dt} = \gamma \ ein + \delta``

#   - ``\gamma=0`` (Linear growing at rate ``\delta``): ``ein(t) = \delta \ t + ein_0``
#   - ``\gamma \neq 0`` (Exponential growing at rate ``\gamma``): ``ein(t) = (ein_0 + \delta / \gamma)e^{\gamma t} - \delta / \gamma``


# ## Reachability settings

using ReachabilityAnalysis, ModelingToolkit
using Plots, Plots.PlotMeasures, LaTeXStrings

LazySets.set_ztol(Float64, 1e-14);

# ``\gamma`` and ``\delta`` control the shape of the input voltage ``ein``

const var = @variables eₒ ein

function opamp_circuit_with_saturation_MT(; X0 = BallInf(zeros(2), 0.0),
                                            R₁ = 2., R₂ = 6., C = 1.e-3,
                                            γ = 100., δ = 100., Es = 2.)
    α = hcat(-1/(R₂ * C))
    β = hcat(-1/(R₁ * C))
    ## transition graph
    automaton = LightAutomaton(2)
    add_transition!(automaton, 1, 2, 1)
    ## modes
    A = [α β; 0. γ]; b = [0.; δ];
    invariant = HalfSpace(ein <= Es, var)
    mode1 = @system(x' = Ax + b, x ∈ invariant)
    A = [α β; 0. 0.];
    mode2 = @system(x' = Ax, x ∈ Universe(2))
    modes = [mode1, mode2]
    ## transition mode1 -> mode2
    guard = Hyperplane([0.; 1.], Es) # ein == Es
    ## transition annotations
    t1 = ConstrainedIdentityMap(2, guard)
    resetmaps = [t1]
    H = HybridSystem(automaton, modes, resetmaps, [AutonomousSwitching()])
    ## initial condition is at the orgin in mode 1
    initial_condition = [(1, X0)]
    return IVP(H, initial_condition)
end

# ## Constant input voltage

# - ``ein = Es``: ``\delta = 0``, ``\gamma = 0``

X0=Hyperrectangle([0., 1.5],[0.0, 0.1]);

const_prob = opamp_circuit_with_saturation_MT(X0=X0, γ = 0., δ = 0., Es = 3.);
const_sol = solve(const_prob, T=0.1, alg=BOX(δ=5e-4));

#-

fig = Plots.plot();
Plots.plot!(fig, const_sol, linecolor=:black, vars=(0, 2),
            color=:blue, alpha=0.4, lw=1.0, label = "",
            xlab=L"t", ylab=L"ein",
            xtick=[0, 0.025, 0.05, 0.075, 0.1], ytick=[0., 0.5, 1.0, 1.5, 2.0],
            xlims=(0., 0.1), ylims=(0., 2.0),
            bottom_margin=6mm, left_margin=2mm, right_margin=4mm, top_margin=3mm,
            legend = :topright
            );
plot!(x->x, x->1.5, 0., 0.1, line=2, color="red", linestyle=:dash, legend=nothing)
fig

#-

fig = Plots.plot();
Plots.plot!(fig, const_sol, linecolor=:black, vars=(0, 1),
            color=:blue, alpha=0.4, lw=1.0, label = "",
            xlab=L"t", ylab=L"e_o",
            xtick=[0, 0.025, 0.05, 0.075, 0.1], ytick=[-5., -4., -3., -2., -1., 0.],
            xlims=(0., 0.1), ylims=(-5., 0.),
            bottom_margin=6mm, left_margin=2mm, right_margin=4mm, top_margin=3mm,
            legend = :topright
            );
plot!(x->x, x->-4.5, 0., 0.1, line=2, color="red", linestyle=:dash, legend=nothing)
fig

# ## Non constant input voltage

# - Linear growing input voltage up to saturation at ``Es``: ``\delta \neq 0``, ``\gamma = 0``
# - Exponential growing input voltage up to saturation at ``Es``: ``\delta \neq 0``, ``\gamma > 0``

X0=Hyperrectangle([0.0, 0.0],[0.0, 0.0]);

lin_prob = opamp_circuit_with_saturation_MT(X0=X0, γ = 0., δ = 100., Es = 3.);
lin_sol = solve(lin_prob, T=0.1, alg=BOX(δ=5e-4));

exp_prob = opamp_circuit_with_saturation_MT(X0=X0, γ = 100., δ = 100., Es = 3.);
exp_sol = solve(exp_prob, T=0.1, alg=BOX(δ=5e-4));

#-

fig = Plots.plot();
Plots.plot!(fig, lin_sol, linecolor=:black, vars=(0, 2),
            color=:blue, alpha=0.4, lw=1.0, label = "Linear");
Plots.plot!(fig, exp_sol, linecolor=:black, vars=(0, 2),
            color=:red, alpha=0.4, lw=1.0, label = "Exponential",
            xlab=L"t", ylab=L"e_o",
            xtick=[0, 0.025, 0.05, 0.075, 0.1], ytick=[0., 1., 2., 3.],
            xlims=(0., 0.1), ylims=(0., 3.5),
            bottom_margin=6mm, left_margin=2mm, right_margin=4mm, top_margin=3mm,
            legend = :bottomright
            );
fig

#-

fig = Plots.plot();
Plots.plot!(fig, lin_sol, linecolor=:black, vars=(0, 1),
            color=:blue, alpha=0.4, lw=1.0, label = "Linear");
Plots.plot!(fig, exp_sol, linecolor=:black, vars=(0, 1),
            color=:red, alpha=0.4, lw=1.0, label = "Exponential",
            xlab=L"t", ylab=L"ein",
            xtick=[0, 0.025, 0.05, 0.075, 0.1], ytick=[-9., -6., -3., 0.],
            xlims=(0., 0.1), ylims=(-10., 0.),
            bottom_margin=6mm, left_margin=2mm, right_margin=4mm, top_margin=3mm,
            legend = :topright
            );
fig

# ## References

# [1] Kluever, C. A. (2015). Dynamic systems: modeling, simulation, and control. John Wiley & Sons.
