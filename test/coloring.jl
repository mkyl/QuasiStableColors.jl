using Test

using Graphs

include("datasets.jl")

"""Returns true if `P` is a valid relative `q`-coloring of G=(`V`,`E`)."""
function verify_q_color(G::AbstractGraph{T}, P::Vector{Vector{T}},
    q::Number) where {T}
    # check disjointness
    if mapreduce(length, +, P) != nv(G)
        @warn "Vertex appears in more than one partition"
        return false
    end

    # check completness
    if reduce(∪, P; init=Set{T}()) != Set(vertices(G))
        @warn "Partition does not cover all nodes in V"
        return false
    end

    max_q = 0

    # check main property
    for P_i in P, P_j in P
        # count number of P_j neighbors
        N = map(x -> length(neighbors(G, x) ∩ P_j), P_i)
        # TODO check inneighbors for digraphs
        max_q = max(max_q, maximum(N) - minimum(N))
    end
    if max_q > q
        @warn "max_q measured at: $max_q"
        return false
    end
    return true
end


@testset "stable coloring correctness" begin
    @testset "one-color chain" begin
        G = SimpleGraph(4)
        add_edge!(G, 1, 1)
        add_edge!(G, 1, 2)
        add_edge!(G, 2, 3)
        add_edge!(G, 3, 4)
        add_edge!(G, 4, 4)

        P = refine_stable(G)
        V = Set(vertices(G))
        @test Set(Set(x) for x in P) == Set([V])
    end

    @testset "two-color chain" begin
        G = SimpleGraph(3)
        add_edge!(G, 1, 2)
        add_edge!(G, 2, 3)

        P = refine_stable(G)
        @test Set(Set(x) for x in P) == Set([Set([1, 3]), Set([2])])
    end

    # see Figure 3, https://doi.org/10.1145/3375395.3387641
    @testset "five colors" begin
        G = SimpleGraph(9)
        add_edge!(G, 1, 5)
        add_edge!(G, 2, 4)
        add_edge!(G, 2, 5)
        add_edge!(G, 2, 8)
        add_edge!(G, 3, 5)
        add_edge!(G, 3, 9)
        add_edge!(G, 6, 9)
        add_edge!(G, 7, 8)
        add_edge!(G, 8, 9)

        P = refine_stable(G)
        P_true = Set([Set([6, 1]), Set([3]), Set([5, 9]), Set([2, 8]), Set([4, 7])])
        @test length(P) == 5
        @test Set(Set(x) for x in P) == P_true
    end
end

@testset "abs approx stable coloring correctness" begin
    @testset "one-color chain" begin
        G = SimpleGraph(4)
        add_edge!(G, 1, 1)
        add_edge!(G, 1, 2)
        add_edge!(G, 2, 3)
        add_edge!(G, 3, 4)
        add_edge!(G, 4, 4)

        P = q_color(G, q=2.0)
        V = Set(vertices(G))
        @test Set(Set(x) for x in P) == Set([V])
    end

    @testset "two-color chain" begin
        G = SimpleGraph(3)
        add_edge!(G, 1, 2)
        add_edge!(G, 2, 3)

        P = q_color(G, q=3.0)
        V = Set(vertices(G))
        @test Set(Set(x) for x in P) == Set([V])
    end

    @testset "OpenFlights" begin
        G = datasets.openflight()
        c = 4.0
        # to test data structure resizing logic
        P = q_color(G, q=c)
        @test verify_q_color(G, P, c)
    end

    @test_skip @testset "DBLP co-authorship" begin
        G = datasets.dblp()
        c = 8.0
        P = q_color(G, q=c)
        @test verify_q_color(G, P, c)
    end
end

@testset "directed graph correctness" begin
    @testset "three-color linked list" begin
        G = SimpleDiGraph(4)
        add_edge!(G, 1, 2)
        add_edge!(G, 2, 3)
        add_edge!(G, 3, 4)

        P = q_color(G, q=0.0)
        P_true = Set([Set([1,]), Set([2, 3]), Set([4,])])
        @test Set(Set(x) for x in P) == P_true
    end

    @testset "one-color linked list" begin
        G = SimpleDiGraph(4)
        add_edge!(G, 1, 2)
        add_edge!(G, 2, 3)
        add_edge!(G, 3, 4)

        P = q_color(G, q=2.0)
        V = Set(vertices(G))
        @test Set(Set(x) for x in P) == Set([V])
    end

    @testset "simple graph" begin
        edges = [
            Edge(1, 3),
            Edge(1, 4),
            Edge(2, 1),
            Edge(4, 1),
        ]
        G = SimpleDiGraphFromIterator(edges)

        P = q_color(G, q=2.0)
        V = Set(vertices(G))
        @test Set(Set(x) for x in P) == Set([V])
    end
end
