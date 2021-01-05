# =========================================================================================================================
# Theoretical error bounds for Carleman linearization
#
# References:
#
# - [1] Forets, Marcelo, and Amaury Pouly. "Explicit error bounds for carleman linearization."
#       arXiv preprint arXiv:1711.02552 (2017).
#
# - [2] Liu, J. P., Kolden, H. Ø., Krovi, H. K., Loureiro, N. F., Trivisa, K., & Childs, A. M. (2020).
#       Efficient quantum algorithm for dissipative nonlinear differential equations. arXiv preprint arXiv:2011.03185.
# =========================================================================================================================

# --- Error bounds using a priori estimate from [1] ---

# See Theorem 4.2 in [1]. This is a bound based on an a priori estimate
# of the norm of the exact solution x(t)
function error_bound_apriori(α, F₁, F₂; N, p=Inf)
    nF₂ = opnorm(F₂, p)
    μF₁ = logarithmic_norm(F₁, p)

    β = α * nF₂ / μF₁
    ε = t -> α * β^N * (exp(μF₁ * t) - 1)^N
    return ε
end

# See Theorem 4.2 in [1]
function convergence_radius_apriori(α, F₁, F₂; N, p=Inf)
    nF₂ = opnorm(F₂, p)
    μF₁ = logarithmic_norm(F₁, p)

    if μF₁ < 0
        return Inf
    end
    β = α * F₂ / μF₁
    T = (1/μF₁) * log(1 + 1/β)
    return T
end

# --- Error bounds using power series method from [1] ---

# See Theorem 4.3 in [1], which uses the power series method
function error_bound_pseries(x₀, F₁, F₂; N, p=Inf)
    nx₀ = norm(x₀, p)
    nF₁ = opnorm(F₁, p)
    nF₂ = opnorm(F₂, p)

    β₀ = nx₀ * nF₂ / nF₁
    ε = t -> nx₀ * exp(nF₁ * t) / (1 - β₀ * (exp(nF₁ * t) - 1)) * (β₀ * (exp(nF₁ * t) - 1))^N
    return ε
end

# See Theorem 4.3 in [1]
function error_bound_apriori(α, F₁, F₂; N, p=Inf)
    nF₂ = opnorm(F₂, p)
    μF₁ = logarithmic_norm(F₁, p)

    β = α * nF₂ / μF₁
    ε = t -> α * β^N * (exp(μF₁ * t) - 1)^N
    return ε
end

function convergence_radius_pseries(x₀, F₁, F₂; N, p=Inf)
    nx₀ = norm(x₀, p)
    nF₁ = opnorm(F₁, p)
    nF₂ = opnorm(F₂, p)

    β₀ = nx₀ * nF₂ / nF₁
    T = (1/nF₁) * log(1 + 1/β₀)
    return T
end

# --- Error bounds using spectral abscissa from [2] ---

# See Lemma 2 in [2]
function error_bound_specabs(x₀, F₁, F₂; N, p=Inf, check=true)
    R = error_bound_specabs_R(x₀, F₁, F₂; p=p)
    if check
        @assert Re_λ₁ < 1 "expected R < 1, got R = $R; try scaling the ODE"
        @assert R < 1 "expected R < 1, got R = $R; try scaling the ODE"
    end
    _error_bound_specabs(nx₀, R, N, Re_λ₁)
end

# See Lemma 2 in [2]
function error_bound_specabs(nx₀, R, N, Re_λ₁)
    ε = t -> nx₀ * R^N * (1 - exp(Re_λ₁ * t))^N
    return ε
end

# See Definition (...) in [2]
function error_bound_specabs_R(x₀, F₁, F₂; p=Inf)
    nx₀ = norm(x₀, p)

    # compute eigenvalues and sort them by increasing real part
    λ = eigvals(F₁, sortby=real)
    λ₁ = last(λ)
    Re_λ₁ = real(λ₁)
    nF₂ = opnorm(F₂, p)
    R = nx₀ * nF₂ / abs(Re_λ₁)
    return R
end

# See Definition (...) in [2]
function error_bound_specabs_zero_jac_R(x₀, F₂; p=Inf)

end

# See Definition (...) in [2]
function convergence_radius_specabs_zero_jac(x₀, F₁, F₂; p=Inf)

end
