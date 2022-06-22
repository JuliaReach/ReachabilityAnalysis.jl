function post(alg::FLOWSTAR, ivp::IVP{<:AbstractContinuousSystem}, timespan;
              Δt0::TimeInterval=zeroI, kwargs...) where {N}

@requires Flowstar

@unpack remainder_estimation, precondition, order_min, order_max, cutoff, precision = alg

# initial time and final time
t0, T = tstart(timespan), tend(timespan)

# vector field
f! = (islinear(ivp) || isaffine(ivp)) ? inplace_field!(ivp) : VectorField(ivp)

# initial set
ivp_norm = _normalize(ivp)
X0 = initial_state(ivp_norm)

# preallocate output flowpipe
F = Vector{TaylorModelReachSet{N}}()

# ...
end
