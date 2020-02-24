function default_discrete_post(problem::InitialValueProblem{ST, XT}) where {ST<:HybridSystem, XT}
    # TODO
end

#=
function _solve_hybrid(problem, post, kwargs)
    # TODO
end
=#

# TODO hy hybrid system ?
function solve(problem::InitialValueProblem{ST},
               opC::AbstractContinuousPost=default_continuous_post(problem),
               opD::AbstractContinuousPost=default_discrete_post(problem);
               T::Float64) where {ST<:HybridSystem}

    # check arguments

    # solve
    _solve_hybrid(problem, opC, opD; T=T)
end


# check arguments
# use arg check... ? paqute ArgCheck
# por el tema de arg check ver HomotopyContinuation, ver tambien el link
# a uno de los issues de MathematicalSystems.
