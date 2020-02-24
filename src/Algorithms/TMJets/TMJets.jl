export TMJets

using Reexport
@reexport using TaylorIntegration

using IntervalArithmetic: IntervalBox

using StaticArrays: SVector

struct TMJets <: AbstractContinuousPost
    options::TwoLayerOptions

    function TMJets(𝑂::Options)
        𝑂new = validate_and_wrap_options(𝑂, options_TMJets())
        return new(𝑂new)
    end
end

# TODO: add invariant in the loop, see
# https://github.com/JuliaReach/Reachability.jl/pull/595/files
# the branch is mforets/property_TMJetspost

# convenience constructor from pairs of symbols
TMJets(𝑂::Pair{Symbol, <:Any}...) = TMJets(Options(Dict{Symbol, Any}(𝑂)))

# default options (they are added in the function validate_and_wrap_options)
TMJets() = TMJets(Options())

include("init.jl")
include("post.jl")
include("project.jl")
include("reach.jl")
