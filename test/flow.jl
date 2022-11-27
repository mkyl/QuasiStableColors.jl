using QuasiStableColors: Flow

using Test
using Graphs
using GraphsFlows

@testset "lifted maxflow, stable" begin
    @testset "simplest 2-flow" begin
        G = SimpleDiGraph(4)
        add_edge!(G, 1, 2)
        add_edge!(G, 1, 3)
        add_edge!(G, 2, 4)
        add_edge!(G, 3, 4)
        @test Flow.lifted_maxflow(G, 1, 4; q=0.0) == 2
    end

    @testset "simple 3-flow" begin
        G = SimpleDiGraph(5)
        add_edge!(G, 1, 2)
        add_edge!(G, 1, 3)
        add_edge!(G, 1, 4)
        add_edge!(G, 2, 5)
        add_edge!(G, 3, 5)
        add_edge!(G, 4, 5)
        @test Flow.lifted_maxflow(G, 1, 5; q=0.0) == 3
    end

    @test_skip @testset "airports" begin
        G, names = datasets.airports()
        s = names["SFO"]
        t = names["JFK"]
        f_ref, _ = maximum_flow(G, s, t)
        @test lifted_maxflow(G, s, t, q=0.0) == f_ref

        s = names["SEA"]
        t = names["JNB"]
        f_ref, _ = maximum_flow(G, s, t)
        @test Flow.lifted_maxflow(G, s, t; q=0.0) == f_ref
    end
end

@testset "lifted maxflow, absolute" begin
    @test_skip @testset "airports SFO to JFK" begin
        G, names = datasets.airports()
        s = names["SFO"]
        t = names["JFK"]
        f_ref, _ = maximum_flow(G, s, t)
        @test 0.9 <= (lifted_maxflow(G, s, t, 2, abs) / f_ref) <= 1.1

        G, names = datasets.airports()
        s = names["SEA"]
        t = names["JNB"]
        f_ref, _ = maximum_flow(G, s, t)
        @test 0.9 <= (lifted_maxflow(G, s, t, 2, abs) / f_ref) <= 1.15
    end
end
