# default continuous post for linear systems
function post(alg::Union{A20,BFFPSV18,BOX,GLGM06,INT,LGG09,VREP},
              ivp::IVP{<:AbstractContinuousSystem}, tspan;
              Δt0::TimeInterval=zeroT, kwargs...)
    δ = alg.δ

    NSTEPS = get(kwargs, :NSTEPS, compute_nsteps(δ, tspan))

    # normalize system to canonical form
    ivp_norm = _normalize(ivp)

    # homogenize system
    if get(kwargs, :homogenize, false)
        ivp_norm = homogenize(ivp_norm)
    end

    # discretize system
    ivp_discr = discretize(ivp_norm, δ, alg.approx_model)

    # apply discrete post to discretized system
    return post(alg, ivp_discr, NSTEPS; Δt0=Δt0, kwargs...)
end
