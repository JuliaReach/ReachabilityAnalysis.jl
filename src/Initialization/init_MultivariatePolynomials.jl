eval(quote
    using .MultivariatePolynomials: AbstractVariable, AbstractMonomialLike,
                                    exponents, variables, powers
end)

eval(load_kron_multivariate())
