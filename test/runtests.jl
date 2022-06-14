using QuasiStableColors
using Test

using Graphs

@testset "QuasiStableColors.jl" begin
    @testset "Coloring" begin
        include("coloring.jl")
    end
    @testset "Linear Optimization" begin
        include("optimize.jl")
    end
    @testset "Max-flow" begin
        include("flow.jl")
    end
    @testset "Centrality" begin
        include("centrality.jl")
    end
end
