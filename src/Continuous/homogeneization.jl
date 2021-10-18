# ====================================================
# Homogeneization of linear systems
# ====================================================

# no-op
homogenize(ivp::IVP{LCS{N, MT}, ST}) where {N, MT<:AbstractMatrix{N}, ST} = ivp
homogenize(sys::LCS{N, MT}) where {N, MT<:AbstractMatrix{N}} = sys

"""
    homogenize(ivp::IVP{CLCCS{N,MT,IdentityMultiple{N},XT,ConstantInput{SI}},ST}) where {N, MT<:AbstractMatrix{N}, XT<:LazySet{N}, SI<:Singleton{N}, ST<:LazySet{N}}

Transform an inhomogeneous linear initial-value problem into an homogeneous one
by introducing auxiliary state variables.

### Input

- `ivp` -- initial-value problem

### Output

Homogeneous initial-value problem.

### Notes

This function transforms the canonical initial-value problem ``x' = Ax + u``,
``x ∈ X`` with ``u(0) ∈ U = {u}`` (a singleton) into an homogeneous problem
without inputs ``y' = Â * y``, ``y ∈ Y``.
"""
function homogenize(ivp::IVP{CLCCS{N,MT,IdentityMultiple{N},XT,ConstantInput{SI}},ST}) where {N, MT<:AbstractMatrix{N}, XT<:LazySet{N}, SI<:Singleton{N}, ST<:LazySet{N}}
    # homogenized state matrix
    U = inputset(ivp) |> ReachabilityAnalysis.next_set
    A = state_matrix(ivp)
    Â = _homogenize_state_matrix(A, U)

    # homogenized input set
    X0 = initial_state(ivp)
    Y0 = _homogenize_initial_state(X0)

    X = stateset(ivp)
    Y = _homogenize_stateset(X)
    ivph = IVP(CLCS(Â, Y), Y0)
    return ivph
end

"""
    homogenize(sys::SOACS)

Transform an inhomogeneous second order system into an homogeneous one
by introducing auxiliary state variables.

### Input

- `sys` -- second order system

### Output

First-order homogeneous system.

### Notes

This function transforms the second-order system ``Mx'' + Cx' + Kx = b`` into a
first-order, homogeneous one, ``y' = Â * y``. It is assumed that the matrix ``M``
is invertible.
"""
function homogenize(sys::SOACS)
    @unpack M, C, K, b = sys

    n = size(M, 1)
    invM = inv(M)
    Zn = spzeros(n, n)
    In = Matrix(1.0I, n, n)

    A = [Zn           In     ;
        -invM*K       -invM*C]
    b0 = vcat(zeros(n), invM * b)

    m = 2n
    Aext = spzeros(m+1, m+1)
    Aext[1:m, 1:m] .= A
    Aext[1:m, m+1] .= b0

    return @system(x' = Aext*x)
end

function _homogenize_state_matrix(A::AbstractMatrix, U::Singleton)
    u = element(U)
    n = size(A, 1)
    Â = zeros(n+1, n+1)
    Â[1:n, 1:n] .= A
    Â[1:n, n+1] .= u
    return Â
end

function _homogenize_initial_state(X0::Singleton{N}) where {N}
    x0 = element(X0)
    y0 = vcat(x0, one(N))
    Y0 = Singleton(y0)
    return Y0
end

function _homogenize_initial_state(X0::LazySet{N}, dims=1) where {N}
    Y0 = X0 × Singleton(ones(N, dims))
    return Y0
end

function _homogenize_stateset(X::Universe, ndims=1)
    return Universe(dim(X) + ndims)
end

function _homogenize_stateset(X::LazySet, ndims=1)
    return X × Universe(ndims)
end

function _homogenize_stateset(X::HalfSpace{N}, ndims=1) where {N}
    a′ = vcat(X.a, fill(zero(N), ndims))
    return HalfSpace(a′, X.b)
end
