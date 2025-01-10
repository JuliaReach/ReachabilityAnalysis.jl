# continuous post
function post(alg::LGG09, ivp::IVP{<:AbstractContinuousSystem}, tspan;
              Δt0::TimeInterval=zeroI, kwargs...)
    δ = alg.δ

    NSTEPS = get(kwargs, :NSTEPS, compute_nsteps(δ, tspan))

    # normalize system to canonical form
    ivp_norm = _normalize(ivp)

    # homogenize the initial-value problem
    if get(kwargs, :homogenize, false)
        ivp_norm = homogenize(ivp_norm)
    end

    # discretize system
    ivp_discr = discretize(ivp_norm, δ, alg.approx_model)

    return post(alg, ivp_discr, NSTEPS; Δt0=Δt0, kwargs...)
end

# discrete post
function post(alg::LGG09{N,AM,VN,TN}, ivp::IVP{<:AbstractDiscreteSystem}, NSTEPS=nothing;
              Δt0::TimeInterval=zeroI, kwargs...) where {N,AM,VN,TN}
    # dimension check
    template = alg.template
    @assert statedim(ivp) == dim(template) "the problems' dimension $(statedim(ivp)) " *
                                           "doesn't match the dimension of the template directions, $(dim(template))"

    @unpack δ, approx_model, template, static, threaded, vars = alg

    if isnothing(NSTEPS)
        if haskey(kwargs, :NSTEPS)
            NSTEPS = kwargs[:NSTEPS]
        else
            throw(ArgumentError("`NSTEPS` not specified"))
        end
    end

    Φ = state_matrix(ivp)
    Ω₀ = initial_state(ivp)
    X = stateset(ivp)

    if alg.sparse # ad-hoc conversion of Φ to sparse representation
        Φ = sparse(Φ)
    end

    cacheval = Val(alg.cache)

    # true <=> there is no input, i.e. the system is of the form x' = Ax, x ∈ X
    got_homogeneous = !hasinput(ivp)

    # NOTE: option :static currently ignored!

    # preallocate output flowpipe
    SN = SubArray{N,1,Matrix{N},Tuple{Base.Slice{Base.OneTo{Int}},Int},true}
    RT = TemplateReachSet{N,VN,TN,SN}
    F = Vector{RT}(undef, NSTEPS)

    if got_homogeneous
        ρℓ = reach_homog_LGG09!(F, template, Ω₀, Φ, NSTEPS, δ, X, Δt0, cacheval, Val(alg.threaded))
    else
        U = inputset(ivp)
        @assert isa(U, LazySet)
        ρℓ = reach_inhomog_LGG09!(F, template, Ω₀, Φ, NSTEPS, δ, X, U, Δt0, cacheval,
                                  Val(alg.threaded))
    end

    return Flowpipe(F, Dict{Symbol,Any}(:sfmat => ρℓ, :alg_vars => vars))
end
