using ReachabilityAnalysis, SparseArrays, JLD2

examples_dir = normpath(@__DIR__, "..", "..", "..", "examples")
HEAT01_path = joinpath(examples_dir, "Heat3D", "HEAT01.jld2")
