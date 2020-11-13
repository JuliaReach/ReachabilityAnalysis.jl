#
# TODO: fix normalization to use getter functions from MathematicalSystems
#
#export add_dimension,
#       next_set,
#       iscontinuoussystem,
#       ishybridsystem

import Base: *
import LazySets.constrained_dimensions

# common aliases for system's names
const LCS = LinearContinuousSystem
const LDS = LinearDiscreteSystem
const CLCS = ConstrainedLinearContinuousSystem
const CLDS = ConstrainedLinearDiscreteSystem
const CLCCS = ConstrainedLinearControlContinuousSystem
const CLCDS = ConstrainedLinearControlDiscreteSystem
const ACS = AffineContinuousSystem
const ADS = AffineDiscreteSystem
const CACCS = ConstrainedAffineControlContinuousSystem
const CACDS = ConstrainedAffineControlDiscreteSystem
const CACS = ConstrainedAffineContinuousSystem
const CADS = ConstrainedAffineDiscreteSystem
const BBCS = BlackBoxContinuousSystem
const CBBCS = ConstrainedBlackBoxContinuousSystem
const CBBCCS = ConstrainedBlackBoxControlContinuousSystem
const SOLCS = SecondOrderLinearContinuousSystem
const SOACS = SecondOrderAffineContinuousSystem
const SOCLCCS = SecondOrderConstrainedLinearControlContinuousSystem
const SOCACCS = SecondOrderConstrainedAffineControlContinuousSystem
const SecondOrderSystem = Union{SOLCS, SOACS, SOCLCCS, SOCACCS}

# continuous systems that are handled by this library
iscontinuoussystem(T::Type{<:AbstractSystem}) = false
iscontinuoussystem(T::Type{<:LCS}) = true
iscontinuoussystem(T::Type{<:CLCS}) = true
iscontinuoussystem(T::Type{<:CLCCS}) = true
iscontinuoussystem(T::Type{<:ACS}) = true
iscontinuoussystem(T::Type{<:CACCS}) = true
iscontinuoussystem(T::Type{<:CACS}) = true
iscontinuoussystem(T::Type{<:BBCS}) = true
iscontinuoussystem(T::Type{<:CBBCS}) = true
iscontinuoussystem(T::Type{<:CBBCCS}) = true
iscontinuoussystem(T::Type{<:SOLCS}) = true
iscontinuoussystem(T::Type{<:SOACS}) = true
iscontinuoussystem(T::Type{<:SOCLCCS}) = true
iscontinuoussystem(T::Type{<:SOCACCS}) = true

# hybrid systems that are handled by this library
ishybridsystem(T::Type{<:AbstractSystem}) = false
ishybridsystem(T::Type{<:HybridSystem}) = true

# systems that are second order
function is_second_order(::ST) where {ST<:AbstractContinuousSystem}
    ST <: SOLCS || ST <: SOACS || ST <: SOCLCCS || ST <: SOCACCS
end

#export LCS, LDS, CLCS, CLDS, CLCCS, CLCDS, CACCS, CACDS, CACS, CADS, IVP, BBCS,
#       CBBCS, CBBCCS

*(M::AbstractMatrix, input::ConstantInput) =  ConstantInput(M * input.U)

# convenience functions
next_set(inputs::ConstantInput) = collect(nextinput(inputs, 1))[1]
next_set(inputs::AbstractInput, state::Int64) = collect(nextinput(inputs, state))[1]

"""
    add_dimension(A::AbstractMatrix, m=1)

Append one or more zero rows and columns to a matrix.

### Input

- `A` -- matrix
- `m` -- (optional, default: `1`) the number of extra dimensions

### Examples

```jldoctest add_dimension_test
julia> A = [0.4 0.25; 0.46 -0.67]
2×2 Array{Float64,2}:
 0.4    0.25
 0.46  -0.67

julia> add_dimension(A)
3×3 Array{Float64,2}:
 0.4    0.25  0.0
 0.46  -0.67  0.0
 0.0    0.0   0.0
```
To append more than one zero row-column, use the second argument `m`:

```jldoctest add_dimension_test
julia> add_dimension(A, 2)
4×4 Array{Float64,2}:
 0.4    0.25  0.0  0.0
 0.46  -0.67  0.0  0.0
 0.0    0.0   0.0  0.0
 0.0    0.0   0.0  0.0
```
"""
function add_dimension(A::AbstractMatrix, m=1)
    n = size(A, 1)
    return vcat(hcat(A, zeros(n, m)), zeros(m, n+m))
end

"""
    add_dimension(X::LazySet, m=1)::LazySet

Adds an extra dimension to a LazySet through a Cartesian product.

### Input

- `X` -- a lazy set
- `m` -- (optional, default: `1`) the number of extra dimensions

### Examples

```jldoctest add_dimension_set
julia> X = BallInf(ones(9), 0.5);

julia> dim(X)
9

julia> Xext = add_dimension(X);

julia> dim(Xext)
10

julia> X = ZeroSet(4);

julia> dim(add_dimension(X))
5

julia> typeof(X)
ZeroSet{Float64}
```

More than one dimension can be added passing the second argument:

```jldoctest add_dimension_set
julia> Xext = add_dimension(BallInf(zeros(10), 0.1), 4);

julia> dim(Xext)
14
```

### Notes

In the special case that the given set is a zero set, instead of cartesian product
a new zero set with extended dimensions is returned.
"""
function add_dimension(X::LazySet, m=1)::LazySet
    return X * ZeroSet(m)
end

function add_dimension(X::ZeroSet, m=1)::ZeroSet
    return ZeroSet(dim(X)+m)
end

"""
    add_dimension(cs, m=1)

Adds an extra dimension to a continuous system.

### Input

- `cs` -- continuous system
- `m` -- (optional, default: `1`) the number of extra dimensions

### Examples

```jldoctest add_dimension_cont_sys
julia> using MathematicalSystems, SparseArrays

julia> A = sprandn(3, 3, 0.5);

julia> X0 = BallInf(zeros(3), 1.0);

julia> s = InitialValueProblem(LinearContinuousSystem(A), X0);

julia> sext = add_dimension(s);

julia> statedim(sext)
4
```

If there is an input set, it is also extended:

```jldoctest add_dimension_cont_sys
julia> U = ConstantInput(Ball2(ones(3), 0.1));

julia> s = InitialValueProblem(ConstrainedLinearControlContinuousSystem(A, Matrix(1.0I, size(A)), nothing, U), X0);

julia> sext = add_dimension(s);

julia> statedim(sext)
4

julia> dim(next_set(inputset(sext)))
4
```

Extending a system with a varying input set:

If there is an input set, it is also extended:

```jldoctest add_dimension_cont_sys
julia> U = VaryingInput([Ball2(ones(3), 0.1 * i) for i in 1:3]);

julia> s = InitialValueProblem(ConstrainedLinearControlContinuousSystem(A, Matrix(1.0I, size(A)), nothing, U), X0);

julia> sext = add_dimension(s);

julia> statedim(sext)
4

julia> dim(next_set(inputset(sext), 1))
4
```

Extending a varing input set with more than one extra dimension:
1] normalize(::AffineContinuousSystem{Float64,Array{Float64,2},Array{Float64,1}}) at /home/mforets/.julia/dev/ReachabilityAnalysis/src/Continuous/normalization.jl:387
```jldoctest add_dimension_cont_sys
julia> sext = add_dimension(s, 7);

julia> statedim(sext)
10

julia> dim(next_set(inputset(sext), 1))
10
```
"""
function add_dimension(cs, m=1)
    Aext = add_dimension(cs.s.A, m)
    X0ext = add_dimension(cs.x0, m)
    if hasmethod(inputset, Tuple{typeof(cs.s)})
        Uext = map(x -> add_dimension(x, m), inputset(cs))
        s = CLCCS(Aext, Matrix(1.0I, size(Aext)), nothing, Uext)
    else
        s = LCS(Aext)
    end
    return IVP(s, X0ext)
end

# accepted types of non-deterministic inputs (non-canonical form)
const UNCF{N} = Union{<:LazySet{N}, Vector{<:LazySet{N}}, <:AbstractInput} where {N}

# accepted types for the state constraints (non-canonical form)
const XNCF{N} = Union{<:LazySet{N}, Nothing} where {N}

"""
    normalize(system::AbstractSystem)

Transform a mathematical system to a normalized (or canonical) form.

### Input

- `system` -- system; it can be discrete or continuous

### Output

Either the same system if it already conforms to a canonical form, or a new
system otherwise.

### Notes

The normalization procedure consists of transforming a given system type into a
"canonical" format that is used internally. More details are given below.

### Algorithm

The implementation of `normalize` exploits `MathematicalSystems`'s' types, which
carry information about the problem as a type parameter.

Homogeneous ODEs of the form ``x' = Ax, x ∈ \\mathcal{X}`` are canonical if the
associated problem is a `ConstrainedLinearContinuousSystem` and `A` is a matrix.
This type does not handle non-deterministic inputs.

Note that a `LinearContinuousSystem` does not consider constraints on the
state-space (such as an invariant); to specify state constraints, use a
`ConstrainedLinearContinuousSystem`. If the passed system is a `LinearContinuousSystem`
(i.e. no constraints) then the normalization fixes a universal set (`Universe`)
as the constraint set.

The generalization to canonical systems with constraints and possibly time-varying
non-deterministic inputs is considered next. These systems are of the form
``x' = Ax + u, u ∈ \\mathcal{U}, x ∈ \\mathcal{X}``. The system type is
`ConstrainedLinearControlContinuousSystem`, where `A` is a matrix, `X` is a set
and `U` is an input, that is, any concrete subtype of `AbstractInput`.

If `U` is not given as an input, normalization accepts either a `LazySet`, or
a vector of `LazySet`s. In these cases, the sets are wrapped around an appropriate
concrete input type.

If the system does not conform to a canonical form, the implementation tries
to make the transformation; otherwise an error is thrown. In particular, ODEs
of the form ``x' = Ax + Bu`` are mapped into ``x' = Ax + u, u ∈ B\\mathcal{U}``,
where now ``u`` has the same dimensions as ``x``.

The transformations described above are analogous in the discrete case, i.e.
``x_{k+1} = A x_k`` and ``x_{k+1} = Ax_{k} + u_k, u_k ∈ \\mathcal{U}, x_k ∈ \\mathcal{X}``
for the linear and affine cases respectively.
"""
function normalize(system::AbstractSystem)
    throw(ArgumentError("the system type $(typeof(system)) is currently not supported"))
end

# convenience extension to IVPs
function normalize(ivp::InitialValueProblem)
    return IVP(normalize(system(ivp)), initial_state(ivp))
end

# "black-box" like systems are not normalized; algorithms should handle this
normalize(S::Union{BBCS, CBBCS, CBBCCS}) = S

# x' = Ax, in the continuous case
# x+ = Ax, in the discrete case
for (L_S, CL_S) in ((:LCS, :CLCS), (:LDS, :CLDS))
    @eval begin
        function normalize(system::$L_S{N, AN}) where {N, AN<:AbstractMatrix{N}}
            n = statedim(system)
            X = Universe(n)
            return $CL_S(state_matrix(system), X)
        end
    end
end

# x' = Ax, x ∈ X in the continuous case
# x+ = Ax, x ∈ X in the discrete case
for CL_S in (:CLCS, :CLDS)
    @eval begin
        function normalize(system::$CL_S{N, AN, XT}) where {N, AN<:AbstractMatrix{N}, M, XT<:XNCF{M}}
            n = statedim(system)
            X = _wrap_invariant(stateset(system), n)
            return $CL_S(state_matrix(system), X)
        end
    end
end

# x' = Ax + Bu, x ∈ X, u ∈ U in the continuous case
# x+ = Ax + Bu, x ∈ X, u ∈ U in the discrete case
for CLC_S in (:CLCCS, :CLCDS)
    @eval begin
        function normalize(system::$CLC_S{N, AN, BN, XT, UT}) where {N, AN<:AbstractMatrix{N}, NN, BN<:AbstractMatrix{NN}, M, XT<:XNCF{M}, UT<:UNCF{M}}
            n = statedim(system)
            X = _wrap_invariant(stateset(system), n)
            U = _wrap_inputs(inputset(system), input_matrix(system))
            $CLC_S(state_matrix(system), I(n, N), X, U)
        end
    end
end

# x' = Ax + b, x ∈ X, u ∈ U in the continuous case
# x+ = Ax + b, x ∈ X, u ∈ U in the discrete case
for (CA_S, CLC_S) in ((:CACS, :CLCCS), (:CADS, :CLCDS))
    @eval begin
        function normalize(system::$CA_S{N, AN, VN, XT}) where {N, AN<:AbstractMatrix{N}, VN<:AbstractVector{N}, XT<:XNCF{N}}
            n = statedim(system)
            X = _wrap_invariant(stateset(system), n)
            U = _wrap_inputs(affine_term(system))
            $CLC_S(state_matrix(system), I(n, N), X, U)
        end
    end
end

# x' = Ax + b in the continuous case
# x+ = Ax + b in the discrete case
for (A_S, CLC_S) in ((:ACS, :CLCCS), (:ADS, :CLCDS))
    @eval begin
        function normalize(system::$A_S{N, AN, VN}) where {N, AN<:AbstractMatrix{N}, VN<:AbstractVector{N}}
            n = statedim(system)
            X = Universe(n)
            U = _wrap_inputs(affine_term(system))
            $CLC_S(state_matrix(system), I(n, N), X, U)
        end
    end
end

# x' = Ax + Bu + c, x ∈ X, u ∈ U in the continuous case
# x+ = Ax + Bu + c, x ∈ X, u ∈ U in the discrete case
for (CAC_S, CLC_S) in ((:CACCS, :CLCCS), (:CACDS, :CLCDS))
    @eval begin
        function normalize(system::$CAC_S{N, AN, BN, VN, XT, UT}) where {N, AN<:AbstractMatrix{N}, BN<:AbstractMatrix{N}, VN<:AbstractVector{N}, XT<:XNCF{N}, UT<:UNCF{N}}
            n = statedim(system)
            X = _wrap_invariant(stateset(system), n)
            U = _wrap_inputs(inputset(system), input_matrix(system), affine_term(system))
            $CLC_S(state_matrix(system), I(n, N), X, U)
        end
    end
end

# fix type inference
function _normalize(ivp::IVP{LCS{N, IdentityMultiple{N}}, Interval{N, IA.Interval{N}}}) where {N}
    return IVP(CLCS(ivp.s.A, Universe(1)), ivp.x0)
end

# ===========================
# Second order systems
# ===========================

# we left the n-dimensional second order system to a 2n-dimensional first
# order system
# if derivatives_last = true (Default)
# we assum that the first n variables correspond to position and the last n
# to velocities as in x̃ = [x, x']
#
# otherwise, if derivatives_last = false
# we assum that the first n variables correspond to velocities and the last n
# to position as in x̃ = [x', x]
function _second_order_linear_matrix(M::AbstractMatrix{N}, C, K; derivatives_last=true) where {N}
    n = size(M, 1)
    M⁻¹ = inv(Matrix(M))
    Idn = Matrix(one(N)*I, n, n)
    Zn = zeros(N, n, n)

    if derivatives_last
        A = [Zn       Idn   ;
             -M⁻¹*K   -M⁻¹*C]
    else
        A = [-M⁻¹*C   -M⁻¹*K;
             Idn          Zn]
    end
    return A, M⁻¹
end

function normalize(system::SOLCS{N}; derivatives_last=true) where {N}
    n = statedim(system)
    M = mass_matrix(system)
    C = viscosity_matrix(system)
    K = stiffness_matrix(system)
    A, _ = _second_order_linear_matrix(M, C, K; derivatives_last=derivatives_last)
    return normalize(LCS(A))
end

function normalize(system::SOACS{N}; derivatives_last=true) where {N}
    n = statedim(system)
    M = mass_matrix(system)
    C = viscosity_matrix(system)
    K = stiffness_matrix(system)

    A, M⁻¹ = _second_order_linear_matrix(M, C, K; derivatives_last=derivatives_last)

    b = affine_term(system)
    if derivatives_last
        c = vcat(zeros(N, n), M⁻¹*b)
    else
        c = vcat(M⁻¹*b, zeros(N, n))
    end

    return normalize(ACS(A, c))
end

function normalize(system::SOCLCCS{N}; derivatives_last=true) where {N}
    n = statedim(system)
    M = mass_matrix(system)
    C = viscosity_matrix(system)
    K = stiffness_matrix(system)
    B = input_matrix(system)
    X = stateset(system)
    U = inputset(system)

    A, M⁻¹ = _second_order_linear_matrix(M, C, K; derivatives_last=derivatives_last)
    if derivatives_last
        B̃ = vcat(zeros(N, n), M⁻¹*B)
    else
        B̃ = vcat(M⁻¹*B, zeros(N, n))
    end
    return normalize(CLCCS(A, B̃, X, U))
end

function normalize(system::SOCACCS{N}; derivatives_last=true) where {N}
    n = statedim(system)
    M = mass_matrix(system)
    C = viscosity_matrix(system)
    K = stiffness_matrix(system)
    B = input_matrix(system)
    d = affine_term(system)
    X = stateset(system)
    U = inputset(system)

    A, M⁻¹ = _second_order_linear_matrix(M, C, K; derivatives_last=derivatives_last)

    if derivatives_last
        B̃ = vcat(zeros(N, n), M⁻¹*B)
        d̃ = vcat(zeros(N, n), M⁻¹*d)
    else
        B̃ = vcat(M⁻¹*B, zeros(N, n))
        d̃ = vcat(M⁻¹*d, zeros(N, n))
    end
    return normalize(CACCS(A, B̃, X, U, d̃))
end

# ---

@inline isidentity(B::IdentityMultiple) = B.M.λ == oneunit(B.M.λ)

_wrap_invariant(X::LazySet, n::Int) = X
_wrap_invariant(X::Nothing, n::Int) = Universe(n)

_wrap_inputs(U::AbstractInput, B::IdentityMultiple) = isidentity(B) ? U : map(u -> B*u, U)
_wrap_inputs(U::AbstractInput, B::AbstractMatrix) = map(u -> B*u, U)

_wrap_inputs(U::LazySet, B::IdentityMultiple) = isidentity(B) ? ConstantInput(U) : ConstantInput(B*U)
_wrap_inputs(U::LazySet, B::AbstractMatrix) = ConstantInput(B*U)

_wrap_inputs(U::Vector{<:LazySet}, B::IdentityMultiple) = isidentity(B) ? VaryingInput(U) : VaryingInput(map(u -> B*u, U))
_wrap_inputs(U::Vector{<:LazySet}, B::AbstractMatrix) = VaryingInput(map(u -> B*u, U))

_wrap_inputs(U::AbstractInput, B::IdentityMultiple, c::AbstractVector) = isidentity(B) ? map(u -> u ⊕ c, U) : map(u -> B*u ⊕ c, U)
_wrap_inputs(U::AbstractInput, B::AbstractMatrix, c::AbstractVector) = map(u -> B*u ⊕ c, U)

_wrap_inputs(U::LazySet, B::IdentityMultiple, c::AbstractVector) = isidentity(B) ? ConstantInput(U ⊕ c) : ConstantInput(B*U ⊕ c)
_wrap_inputs(U::LazySet, B::AbstractMatrix, c::AbstractVector) = ConstantInput(B*U ⊕ c)

_wrap_inputs(U::Vector{<:LazySet}, B::IdentityMultiple, c::AbstractVector) = isidentity(B) ? VaryingInput(map(u -> u ⊕ c, U)) : VaryingInput(map(u -> B*u ⊕ c, U))
_wrap_inputs(U::Vector{<:LazySet}, B::AbstractMatrix, c::AbstractVector) = VaryingInput(map(u -> B*u ⊕ c, U))

_wrap_inputs(c::AbstractVector) = ConstantInput(Singleton(c))

# TODO refactor to MathematicalSystems.jl ?
function MathematicalSystems.ConstrainedLinearControlContinuousSystem(A::Matrix{<:IntervalArithmetic.Interval}, B::SparseMatrixCSC, X::XT, U::UT) where {XT, UT}
    CLCCS(IntervalMatrix(A), IntervalMatrix(B), X, U)
end
function MathematicalSystems.ConstrainedLinearControlDiscreteSystem(A::IntervalMatrix, B::Matrix, X::XT, U::UT) where {XT, UT}
    CLCDS(A, IntervalMatrix(B), X, U)
end

# ==========================================================
# Shared functionality for linear continuous post operators
# ==========================================================

Base.:*(im::IdentityMultiple, d::Diagonal) = im.M.λ * d  # TODO: add to MathematicalSystems

hasinput(S::AbstractSystem) = inputdim(S) > 0
isconstantinput(::ConstantInput) = true
isconstantinput(::VaryingInput) = false
isconstantinput(::LazySet) = true

# The canonical form is:
#     - If the system doesn't have input, a constrained linear continous system (CLCS)
#       x' = Ax, x ∈ X
#     - If the system has an input, a CLCCS, x' = Ax + u, x ∈ X, u ∈ U
# If the original system is unconstrained, the constraint set X is the universal set.
abstract type AbstractLinearContinuousSystem <: AbstractContinuousSystem end
abstract type AbstractNonlinearContinuousSystem <: AbstractContinuousSystem end

function _normalize(ivp::IVP{<:AbstractContinuousSystem})
    if islinear(ivp) || isaffine(ivp)
        return _normalize(ivp, AbstractLinearContinuousSystem)
    else
        return _normalize(ivp, AbstractNonlinearContinuousSystem)
    end
end

function _normalize(ivp::IVP{<:BBCS}, ::Type{AbstractNonlinearContinuousSystem})
    f = ivp.s.f # TODO add getter function in MathematicalSystems
    n = statedim(ivp)
    X = Universe(n)
    return IVP(CBBCS(f, n, X), initial_state(ivp))
end

function _normalize(ivp::IVP{<:CBBCS}, ::Type{AbstractNonlinearContinuousSystem})
    return ivp
end

function _normalize(ivp::IVP{<:AbstractContinuousSystem}, ::Type{AbstractNonlinearContinuousSystem})
    throw(ArgumentError("can't normalize this nonlinear initial-value problemof type $(typeof(ivp))"))
end

_normalize_initial_state(ivp, X0) = X0
_normalize_initial_state(ivp, X0::AbstractVector) = Singleton(X0)
_normalize_initial_state(ivp, X0::IA.Interval) = convert(Interval, X0)
_normalize_initial_state(ivp, X0::IA.IntervalBox) = convert(Hyperrectangle, X0)
_normalize_initial_state(ivp, X0::Number) = Singleton([X0])
_normalize_initial_state(ivp::IVP{<:SecondOrderSystem}, X0::Tuple{<:AbstractVector, <:AbstractVector}) = Singleton(vcat(X0[1], X0[2]))
_normalize_initial_state(ivp::IVP{<:SecondOrderSystem}, X0::Tuple{<:LazySet, <:LazySet}) = X0[1] × X0[2]

const CanonicalLinearContinuousSystem = Union{CLCS, CLCCS}

function _normalize(ivp::IVP{<:AbstractContinuousSystem}, ::Type{AbstractLinearContinuousSystem})

    # initial states normalization
    X0 = initial_state(ivp)
    X0_norm = _normalize_initial_state(ivp, X0)

    # system's normalization
    S = system(ivp)
    S_norm = normalize(S)

    if S_norm === S && X0_norm === X0
        ivp_norm = ivp
    else
        ivp_norm = IVP(S_norm, X0_norm)
    end
    return ivp_norm
end
