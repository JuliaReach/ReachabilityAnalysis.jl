function __init__()
    @require DynamicPolynomials = "7c1d4256-1411-5781-91ec-d7bc3513ac07" include("init_DynamicPolynomials.jl")
    @require MultivariatePolynomials = "102ac46a-7ee4-5c85-9060-abc95bfdeaa3" include("init_MultivariatePolynomials.jl")
    @require Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7" nothing
end
