using Documenter, ReachabilityAnalysis

DocMeta.setdocmeta!(ReachabilityAnalysis, :DocTestSetup,
                   :(using ReachabilityAnalysis); recursive=true)

# pass --fast as an argument to skip rebuilding the examples and running doctests
const _FAST = findfirst(isequal("--fast"), ARGS) !== nothing

# ========================
# Generate examples
# ========================
if !_FAST
    include("generate.jl")
end

# ========================
# Tutorials
# ========================

SET_REPRESENTATIONS = ["Set representations" => [
    "Introduction"                     => "tutorials/set_representations/introduction.md",
    "Hyperrectangular sets"            => "tutorials/set_representations/hyperrectangles.md",
    "Zonotopic sets"                   => "tutorials/set_representations/zonotopes.md",
    "Polyhedral computations"          => "tutorials/set_representations/polyhedral_computations.md",
    "Lazy set representations"         => "tutorials/set_representations/lazy_sets.md",
    "Computing with support functions" => "tutorials/set_representations/support_functions.md",
    #"Other set operations"             => "tutorials/set_representations/other_set_operations.md",
    "LazySets type hierarchy"          => "tutorials/set_representations/lazysets_hierarchy.md",
    "Numerical tolerance"              => "tutorials/set_representations/tolerance.md",
    "Metric notions"                   => "tutorials/set_representations/distances.md"]]

LINEAR_METHODS = ["Linear methods" => [
    "Quickstart"                       => "tutorials/linear_methods/introduction.md",
    "Discrete time reachability"       => "tutorials/linear_methods/discrete_time.md",
    "Dense time reachability"          => "tutorials/linear_methods/dense_time.md",
#    "Propagating zonotopes"            => "tutorials/linear_methods/zonotopes_setprop.md",
#    "Propagating hyperrectangles"      => "tutorials/linear_methods/hyperrectangles_setprop.md",
#    "Propagating support functions"    => "tutorials/linear_methods/supfunc_setprop.md",
    #"Helicopter model"                 => "models/Helicopter.md",
    "Helicopter model"                 => "tutorials/linear_methods/helicopter.md",
    "Structural model"                 => "models/ISS.md"]]

UNCERTAIN_INPUTS = ["Modeling uncertain inputs" => [
    "Introduction"                     => "tutorials/uncertain_inputs/introduction.md",
    "Building model"                   => "models/Building.md",
    "Operational amplifier"            => "models/OpAmp.md",
    "Transmision line"                 => "models/TransmissionLine.md",]]

TAYLOR_METHODS = ["Taylor methods"     => [
    "Introduction"                     => "tutorials/taylor_methods/introduction.md",
    "Taylor model reach-sets"          => "tutorials/taylor_methods/taylor_model_reachsets.md",
    "Domain splitting"                 => "tutorials/taylor_methods/domain_splitting.md",
    "Common gotchas"                   => "tutorials/taylor_methods/gotchas.md",
    "Lotka-Volterra"                   => "models/LotkaVolterra.md",
    "Van der Pol oscillator"           => "models/VanDerPol.md",
    "Lorenz system"                    => "models/Lorenz.md"]]

UNCERTAIN_PARAMETERS = ["Modeling uncertain parameters" => [
    "Introduction"                     => "tutorials/parametric_reachability/introduction.md",
    ]]

HYBRID_SYSTEMS = ["Hybrid systems" => [
    "Introduction"                     => "tutorials/hybrid_systems/introduction.md",
    "Square wave oscillator"           => "models/SquareWaveOscillator.md",
    "Thermostat model"                 => "tutorials/hybrid_systems/thermostat.md"]]

CLOCKED_SYSTEMS = ["Clocked systems"  => [
    "Introduction"                     => "tutorials/clocked_systems/introduction.md",
    "Platoon"                          => "models/Platoon.md",
    "Electro-mechanic break"           => "models/EMBrake.md"]]

BACKWARDS_REACHABILITY = ["Backwards reachability" => [
    "Introduction"                     => "tutorials/backwards_reachability/introduction.md",
    "Projectile motion"                => "models/Projectile.md"]]

LINEAR_PDE = ["Linear PDEs" => [
    "Introduction"                     => "tutorials/linear_pde/introduction.md",
    "Clamped beam"                     => "tutorials/linear_pde/clamped.md",
    "Heat transfer"                    => "tutorials/linear_pde/heat_transfer.md",
    "Concrete heat of hydration"       => "tutorials/linear_pde/concrete.md"]]
#   "Clamped beam"                     => "examples/Clamped.md",
#   "Heat transfer"                    => "examples/Heat1D.md",
#   "Concrete heat of hydration"       => "examples/Concrete.md"

LINEAR_SOLVERS = ["ASB07"    => "man/algorithms/ASB07.md",
                  "BFFPSV18" => "man/algorithms/BFFPSV18.md",
                  "BOX"      => "man/algorithms/BOX.md",
                  "GLGM06"   => "man/algorithms/GLGM06.md",
                  "INT"      => "man/algorithms/INT.md",
                  "LGG09"    => "man/algorithms/LGG09.md",
                  "ORBIT"    => "man/algorithms/ORBIT.md",
                  "VREP"     => "man/algorithms/VREP.md"]

LINEAR_SOLVERS_API = ["ASB07"    => "lib/algorithms/ASB07.md",
                      "BFFPSV18" => "lib/algorithms/BFFPSV18.md",
                      "BOX"      => "lib/algorithms/BOX.md",
                      "GLGM06"   => "lib/algorithms/GLGM06.md",
                      "INT"      => "lib/algorithms/INT.md",
                      "LGG09"    => "lib/algorithms/LGG09.md",
                      "ORBIT"    => "lib/algorithms/ORBIT.md",
                      "VREP"     => "lib/algorithms/VREP.md"]

NONLINEAR_SOLVERS = ["CARLIN"   => "man/algorithms/CARLIN.md",
                     "QINT"     => "man/algorithms/QINT.md",
                     "TMJets"   => "man/algorithms/TMJets.md"]

NONLINEAR_SOLVERS_API = ["CARLIN"   => "lib/algorithms/CARLIN.md",
                         "QINT"     => "lib/algorithms/QINT.md",
                         "TMJets"   => "lib/algorithms/TMJets.md"]

# ========================
# Docs contents
# ========================

makedocs(
    format = Documenter.HTML(prettyurls = haskey(ENV, "GITHUB_ACTIONS"),
                             collapselevel = 1,
                             assets = ["assets/aligned.css"]),
    sitename = "ReachabilityAnalysis.jl",
    doctest = false,
    strict = false,
    pages = [
        "Overview" => "index.md",

        "Tutorials"    => vcat(SET_REPRESENTATIONS,
                          LINEAR_METHODS,
                          UNCERTAIN_INPUTS,
                          TAYLOR_METHODS,
                          UNCERTAIN_PARAMETERS,
                          LINEAR_PDE,
                          HYBRID_SYSTEMS,
                          CLOCKED_SYSTEMS),

        "Manual"       => ["Systems"           => "man/systems.md",
                           "Reach-sets"        => "man/reachsets.md",
                           "Flowpipes"         => "man/flowpipes.md",
                           "Linear solvers"    => LINEAR_SOLVERS,
                           "Nonlinear solvers" => NONLINEAR_SOLVERS,
                           "Solutions"         => "man/solutions.md",
                           "Invariants"        => "man/invariants.md",
                           "Visualization"     => "man/visualization.md",
                           "Projections"       => "man/projections.md",
                           "Clustering"        => "man/clustering.md"],

        "Benchmarks" => ["Repeatability evaluations"   => "man/benchmarks/repeatability.md",
                         "Model library"               => "man/benchmarks/model_library.md",
                         "Comparison with other tools" => "man/benchmarks/comparison.md",
                         "Benchmark repository"        => "man/benchmarks/benchmarks.md",
                         "Filtered oscillator"         => "man/benchmarks/filtered_oscillator.md"],

        "Further examples" => ["Overview"               => "man/examples_overview.md",
                               "Duffing oscillator"     => "models/DuffingOscillator.md",
                               "Laub-Loomis"            => "models/LaubLoomis.md",
                               "Production-Destruction" => "models/ProductionDestruction.md",
                               "Brusselator"            => "models/Brusselator.md",
                               "Quadrotor"              => "models/Quadrotor.md",
                               "Epidemic model"         => "models/SEIR.md",
                               "Spacecraft"             => "models/Spacecraft.md"],

        "API Reference" => ["Systems"                       => "lib/systems.md",
                            "Reach-sets"                    => "lib/reachsets.md",
                            "Flowpipes"                     => "lib/flowpipes.md",
                            "Linear solvers"                => LINEAR_SOLVERS_API,
                            "Nonlinear solvers"             => NONLINEAR_SOLVERS_API,
                            "Solutions"                     => "lib/solutions.md",
                            "Discretization"                => "lib/discretize.md",
                            "Projections"                   => "lib/projections.md",
                            "Clustering"                    => "lib/clustering.md",
                            "Further set operations"        => "lib/operations.md",
                            "Distributed computations"      => "lib/distributed.md",
                            "Internal functions and macros" => "lib/internals.md"],

        "Frequently Asked Questions (FAQ)" => "man/faq.md",

        "Bibliography" => "references.md",

        "How to contribute" => "about.md"
    ]
)

# Deploy built documentation from Travis.
deploydocs(
    repo = "github.com/JuliaReach/ReachabilityAnalysis.jl.git",
    push_preview = true,
)
