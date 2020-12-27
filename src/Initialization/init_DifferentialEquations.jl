using .OrdinaryDiffEq

# extend the solve API for initial-value problems
function OrdinaryDiffEq.solve(ivp::InitialValueProblem, args...; kwargs...)
    tspan = _get_tspan(args...; kwargs...)
    tspan = (inf(tspan), sup(tspan))
    inplace = haskey(kwargs, :inplace) ? kwargs[:inplace] : true
    if inplace
        field = inplace_field!(ivp)
    else
        field = outofplace_field(ivp)
    end
    odeprob = ODEProblem(field, initial_state(ivp), tspan)
    return OrdinaryDiffEq.solve(odeprob, args...; kwargs...)
end
