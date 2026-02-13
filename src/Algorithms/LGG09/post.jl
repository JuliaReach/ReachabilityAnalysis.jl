# default continuous post

# discrete post
function post(alg::LGG09{N,AM,VN,TN}, ivp::IVP{<:AbstractDiscreteSystem}, NSTEPS=nothing;
              Δt0::TimeInterval=zeroT, kwargs...) where {N,AM,VN,TN}
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
    RT = rsetrep(alg)
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
