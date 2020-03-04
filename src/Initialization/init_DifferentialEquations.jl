using .DifferentialEquations

function outofplace_field(ivp::InitialValueProblem)
    # function closure over the inital-value problem
    f = function f_outofplace(x, p, t)
             VectorField(ivp)(x)
         end
    return f
end

function inplace_field!(ivp::InitialValueProblem)
    # function closure over the inital-value problem
    f! = function f_inplace!(dx, x, p, t)
             dx = VectorField(ivp)(x)
         end
    return f_inplace!
end

# extend the solve API for initial-value problems
function DifferentialEquations.solve(ivp::InitialValueProblem, args...; kwargs...)
    tspan = _get_tspan(args...; kwargs...)
    tspan = (inf(tspan), sup(tspan))
    inplace = haskey(kwargs, :inplace) ? kwargs[:inplace] : false
    if inplace
        @warn("the inplace option is not working properly")
        field = inplace_field!(ivp)
    else
        field = outofplace_field(ivp)
    end
    odeprob = ODEProblem(field, initial_state(ivp), tspan)
    return DifferentialEquations.solve(odeprob, args...; kwargs...)
end
