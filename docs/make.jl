ENV["GKSwstype"] = "100"  # set 'GR environment' to 'no output' (for Travis CI)
using Documenter, ReachabilityAnalysis

DocMeta.setdocmeta!(ReachabilityAnalysis, :DocTestSetup,
                   :(using ReachabilityAnalysis); recursive=true)

# generate Literate documentation
include("generate.jl")

makedocs(
    format = Documenter.HTML(prettyurls = haskey(ENV, "GITHUB_ACTIONS"),  # disable for local builds
                             collapselevel = 1,
                             assets = ["assets/juliareach-light.css"]),
    sitename = "ReachabilityAnalysis.jl",
    doctest = false,
    strict = false,
    pages = [
        "Home" => "index.md",
        "Manual" => Any["Basics" => "man/basics.md",
                        "Linear ODEs" => "man/linear.md",
                        "Nonlinear ODEs" => "man/nonlinear.md",
                        "Semidiscrete PDEs" => "man/pde.md",
                        "Hybrid systems" => "man/hybrid.md",
                        "Clocked systems" => "man/clocked.md",
                        "Exploiting structure" => "man/structure.md",
                        "Parametric reachability" => "man/parametric.md",
                        "Backwards reachability" => "man/backwards.md",
                        "Benchmarks" => "man/benchmarks.md",
                        "Model library" => "man/library.md",
                        "Frequently Asked Questions" => "man/faq.md"],
        "Examples" => Any["Overview" => "man/examples_overview.md",
                          "Operational amplifier" => "models/OpAmp.md",
                          "Square wave oscillator" => "models/SquareWaveOscillator.md",
                          "Van der Pol oscillator" => "models/VanDerPol.md",
                          "Duffing oscillator" => "models/DuffingOscillator.md",
                          "Transmision line" => "models/TransmissionLine.md",
                          "Laub-Loomis" => "models/LaubLoomis.md",
                          "Building" => "models/Building.md",
                          "Production-Destruction" => "models/ProductionDestruction.md",
                          "Lotka-Volterra" => "models/LotkaVolterra.md",
                          "Brusselator" => "models/Brusselator.md",
                          "ISS" => "models/ISS.md",
                          "Lorenz system" => "models/Lorenz.md",
                          "Platoon" => "models/Platoon.md",
                          "Quadrotor" => "models/Quadrotor.md",
                          "SEIR model" => "models/SEIR.md",
                          "Spacecraft" => "models/Spacecraft.md"],
                          # "Electromechanic break" => "man/applications/EMBrake.md"
        "Algorithms" => Any["ASB07" => "lib/algorithms/ASB07.md",
                            "BFFPSV18" => "lib/algorithms/BFFPSV18.md",
                            "BOX" => "lib/algorithms/BOX.md",
                            "CARLIN" => "lib/algorithms/CARLIN.md",
                            "GLGM06" => "lib/algorithms/GLGM06.md",
                            "INT" => "lib/algorithms/INT.md",
                            "LGG09" => "lib/algorithms/LGG09.md",
                            "ORBIT" => "lib/algorithms/ORBIT.md",
                            "QINT" => "lib/algorithms/QINT.md",
                            "TMJets" => "lib/algorithms/TMJets.md",
                            "VREP" => "lib/algorithms/VREP.md"],
                            #
        "API Reference" => Any["Reach-sets" => "lib/reachsets.md",
                               "Flowpipes" => "lib/flowpipes.md",
                               "Solutions" => "lib/solutions.md",
                               "Systems" => "lib/systems.md",
                               "Discretization" => "lib/discretize.md",
                               "Projections" => "lib/projections.md",
                               "Clustering" => "lib/clustering.md",
                               "Further set operations" => "lib/operations.md",
                               "Distributed computations" => "lib/distributed.md",
                               "Internal functions and macros" => "lib/internals.md"],
        "References" => "references.md",
        "About" => "about.md"
    ]
)

# Deploy built documentation from Travis.
deploydocs(
    repo = "github.com/JuliaReach/ReachabilityAnalysis.jl.git",
    push_preview = true,
)
