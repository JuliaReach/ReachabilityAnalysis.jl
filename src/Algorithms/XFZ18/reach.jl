function load_XFZ18_reach()
    return quote
        const ∂ = differentiate

        # build the SOS problem
        function build_sos(𝒮, opt)
            T = opt[:T]

            # scale dynamics
            T = opt[:T]
            f = T * 𝒮.s.p

            # define polynomial symbolic variables
            vars = variables(𝒮.s) # TODO: variables(𝒮)
            @polyvar t

            k = opt[:relaxation_degree]

            # monomial vector up to order k
            # 0 <= sum_i alpha_i <= k, if alpha_i is the exponent of x_i
            X = monomials(vars, 0:k)
            XT = monomials([vars; t], 0:k)

            # create a SOS JuMP model to solve with Mosek
            solver = opt[:solver]
            model = SOSModel(with_optimizer(solver))

            # add unknown Φ to the model
            @variable(model, Φ, Poly(XT))

            # jacobian
            ∂t = α -> ∂(α, t)
            ∂xf = α -> ∂(α, x₁) * f[1] + ∂(α, x₂) * f[2]
            LΦ = ∂t(Φ) + ∂xf(Φ)

            # Φ(x, t) at time 0
            Φ₀ = subs(Φ, t => 0.0)

            # scalar variable
            @variable(model, ϵ)

            dom1 = @set t * (T - t) >= 0 && g >= 0
            dom2 = @set g >= 0
            @constraint(model, ϵ >= 0.0)
            @constraint(model, LΦ ∈ SOSCone(), domain=dom1)
            @constraint(model, ϵ - LΦ ∈ SOSCone(), domain=dom1)
            @constraint(model, Φ₀ - V₀ ∈ SOSCone(), domain=dom2)
            @constraint(model, ϵ + V₀ - Φ₀ ∈ SOSCone(), domain=dom2)

            @objective(model, Min, ϵ)
            return model
        end

        # solve model, check feasibility and return polynomials
        function solve_sos(model; verbose=true)
            optimize!(model)

            if verbose
                println("Relaxation order : k = $k")
                println("JuMP.termination_status(model) = ", JuMP.termination_status(model))
                println("JuMP.primal_status(model) = ", JuMP.primal_status(model))
                println("JuMP.dual_status(model) = ", JuMP.dual_status(model))
                println("JuMP.objective_bound(model) = ", JuMP.objective_bound(model))
                println("JuMP.objective_value(model) = ", JuMP.objective_value(model))
            end

            # TODO: error if it fails?
        end

        function extract_approximations(model, 𝑂)
            # time horizon TODO : check consistency w/rescaling
            T = 𝑂[:T]

            # Recovering the solution:
            ϵopt = JuMP.objective_value(model)

            # Punder <= 0  TODO: @set Punder <= 0 ?
            Punder = subs(JuMP.value(model[:Φ]), t => T)

            # Pover <= 0   TODO: @set Pover <= 0 ?
            Pover = subs(JuMP.value(model[:Φ]), t => T) - ϵopt * (T + 1)

            return (ϵopt, Punder, Pover)
        end
    end
end  # quote / load_XFZ18_reach()

function post(::XFZ18, 𝒮::PolynomialContinuousSystem, 𝑂)
    # construct sum-of-squares problem
    model = build_sos(𝒮, opt)

    # solve the sum-of-squares optimization
    solve_sos(model)

    # extract under and over approximations
    ϵopt, Punder, Pover = extract_approximations(model, 𝑂)

    # returns the polynomial under and overapproximations of the reach set
    # for any t ∈ [0, T]
    return (Punder, Pover)
end
