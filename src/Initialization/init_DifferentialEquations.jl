using .DifferentialEquations

function outofplace_field(ivp::LinearContinuousSystem)
    # function closure over the inital-value problem
    f = function outofplace_field(x, p, t)
             VectorField(ivp)(x)
         end
    return f
end

function inplace_field!(ivp::LinearContinuousSystem)
    # function closure over the inital-value problem
    f! = function inplace_field!(dx, x, p, t)
             dx = VectorField(ivp)(x)
         end
    return f!
end

# extend the solve API for initial-value problems
function DifferentialEquations.solve(ivp::InitialValueProblem, args...; kwargs...)
    tspan = _get_tspan(args...; kwargs...)
    odeprob = ODEProblem(inplace_field!(ivp), initial_state(ivp), tspan) # choose between inplace or outoflace
    return DifferentialEquations.solve(odeprob, args...; kwargs...)
end
