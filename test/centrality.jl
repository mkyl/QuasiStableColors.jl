using QuasiStableColors: Centrality

using Test
using Graphs

@testset "betweenness centrality" begin
    @testset "stable-coloring 4 node" begin
        E = [Edge(1, 2), Edge(1, 3), Edge(2, 4), Edge(3, 4), Edge(4, 5)]
        G = SimpleGraphFromIterator(E)
        C = betweenness_centrality(G, normalize=false)
        C₀ = Centrality.approx_betweenness_centrality(G, eps=0.0)
        @test C == C₀
    end

    @testset "stable-coloring 4 node, self-loop" begin
        E = [Edge(1, 2), Edge(1, 3), Edge(2, 3), Edge(2, 4), Edge(3, 4), Edge(4, 5)]
        G = SimpleGraphFromIterator(E)
        C = betweenness_centrality(G, normalize=false)
        C₀ = Centrality.approx_betweenness_centrality(G, eps=0.0)
        @test C == C₀
    end

    # see Figure 3, https://doi.org/10.1145/3375395.3387641
    @test_broken @testset "stable-coloring five colors" begin
        E = [Edge(1, 5), Edge(2, 4), Edge(2, 5), Edge(2, 8), Edge(3, 5), Edge(3, 9),
            Edge(6, 9), Edge(7, 8), Edge(8, 9)]
        G = SimpleGraphFromIterator(E)

        C₀ = Centrality.approx_betweenness_centrality(G, eps=0.0)
        C = betweenness_centrality(G, normalize=false)

        @test_broken C == C₀
    end
end
