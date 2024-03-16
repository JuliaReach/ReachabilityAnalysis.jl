function _template(; n, N)
    dirs = Vector{SingleEntryVector{Float64}}()
    d = sum(n^i for i in 1:N)

    for i in 1:n
        x = SingleEntryVector(i, d, 1.0)
        push!(dirs, x)
    end
    for i in 1:n
        x = SingleEntryVector(i, d, -1.0)
        push!(dirs, x)
    end
    return CustomDirections(dirs)
end

# TODO check if it can be removed
function _project(sol, vars)
    return πsol_1n = Flowpipe([ReachSet(set(project(R, vars)), tspan(R)) for R in sol])
end

struct CanonicalQuadraticForm{T,MT<:AbstractMatrix{T}} <: AbstractContinuousSystem
    F1::MT
    F2::MT
end
statedim(f::CanonicalQuadraticForm) = size(f.F1, 1)
canonical_form(s::CanonicalQuadraticForm) = s.F1, s.F2

# TODO generalize to AbstractContinuousSystem using vector_field
function canonical_form(s::BlackBoxContinuousSystem)
    @requires Symbolics
    return f = s.f
    # differentiate
end

function post(alg::CARLIN, ivp::IVP{<:AbstractContinuousSystem}, tspan; Δt0=interval(0), kwargs...)
    @unpack N, compress, δ, bloat, resets = alg

    # for now we assume there are no resets
    if haskey(kwargs, :NSTEPS)
        NSTEPS = kwargs[:NSTEPS]
        T = NSTEPS * δ
    else
        # get time horizon from the time span imposing that it is of the form (0, T)
        T = _get_T(tspan; check_zero=true, check_positive=true)
        NSTEPS = ceil(Int, T / δ)
    end

    # extract initial states
    X0 = initial_state(ivp)
    # compute canonical quadratic form from the given system
    F1, F2 = canonical_form(system(ivp))
    n = statedim(ivp)
    dirs = _template(; n=n, N=N)
    alg = LGG09(; δ=δ, template=dirs, approx_model=Forward())

    return reach_CARLIN_alg(X0, F1, F2; alg, resets, N, T, Δt0, bloat, compress)
end
