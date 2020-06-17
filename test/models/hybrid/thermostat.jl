# ============================================================================
# Thermostat model
# See https://fenix.tecnico.ulisboa.pt/downloadFile/3779579688470/lygeros.pdf,
# Section 1.3.4.
# ============================================================================

using Reachability, HybridSystems, MathematicalSystems, LazySets, LinearAlgebra
using LazySets: HalfSpace  # resolve name-space conflicts with Polyhedra

function thermostat()
    c_a = 0.1

    # automaton structure
    automaton = LightAutomaton(2)

    # mode on
    A = hcat(-c_a)
    B = hcat(30.)
    U = Singleton([c_a])
    inv = HalfSpace([1.0], 22.0)  # x ≤ 22
    # @vars x v
    # @set x ≤ 22
    m_on = ConstrainedLinearControlContinuousSystem(A, B, inv, ConstantInput(U))

    # mode off
    A = hcat(-c_a)
    B = hcat(0.0)
    U = Singleton([0.0])
    inv = HalfSpace([-1.0], -18.0)  # x ≥ 18
    m_off = ConstrainedLinearControlContinuousSystem(A, B, inv, ConstantInput(U))

    # modes
    modes = [m_on, m_off]

    # transition from on to off
    add_transition!(automaton, 1, 2, 1)
    guard = HalfSpace([-1.0], -21.0)  # x ≥ 21
    t_on2off = ConstrainedIdentityMap(2, guard)

    # transition from off to on
    add_transition!(automaton, 2, 1, 2)
    guard = HalfSpace([1.0], 19.0)  # x ≤ 19
    t_off2on = ConstrainedIdentityMap(2, guard)

    # transition annotations
    resetmaps = [t_on2off, t_off2on]

    # switching
    switchings = [AutonomousSwitching()]

    ℋ = HybridSystem(automaton, modes, resetmaps, switchings)

    # initial condition in mode off
    X0 = Singleton([18.0])
    initial_condition = [(2, X0)]

    system = InitialValueProblem(ℋ, initial_condition)

    options = Options(:mode=>"reach", :T=>5.0,
                      :max_jumps=>1, :verbosity=>1)

    return (system, options)
end
