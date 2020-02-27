#=
function _default_dpost(problem::InitialValueProblem{ST, XT}) where {ST<:HybridSystem, XT}
    # TODO
end

function _solve_hybrid(problem, post, kwargs)
    # TODO
end

# TODO: hybrid system ..
function solve(problem::InitialValueProblem{ST},
               opC::AbstractContinuousPost=default_continuous_post(problem),
               opD::AbstractContinuousPost=default_discrete_post(problem);
               T::Float64) where {ST<:HybridSystem}

    # check arguments

    # solve
    _solve_hybrid(problem, opC, opD; T=T)
end

=#

# check arguments
# use arg check... ? paqute ArgCheck
# por el tema de arg check ver HomotopyContinuation, ver tambien el link
# a uno de los issues de MathematicalSystems.

# ==================
# Argument handling
# ==================

"""
    distribute_initial_set(system::InitialValueProblem{<:HybridSystem, <:LazySet)

Distribute the set of initial states to each mode of a hybrid system.

### Input

- `system` -- an initial value problem wrapping a mathematical system (hybrid)
              and a set of initial states

### Output

A new initial value problem with the same hybrid system but where the set of initial
states is the list of tuples `(state, X0)`, for each state in the hybrid system.
"""
function distribute_initial_set(system::InitialValueProblem{<:HybridSystem, <:LazySet})
    HS, X0 = system.s, system.x0
    initial_states = [(loc, X0) for loc in states(HS)]
    return InitialValueProblem(HS, initial_states)
end
