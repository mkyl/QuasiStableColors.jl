using Graphs
using StatsBase

include("main.jl")

function lifted_centrality(G; args...)
    @assert !is_directed(G)
    G_q, W, _, P = quotient_graph(G; args...)
    C_q = quotient_centrality(G_q, W, P, G)
    C = zeros(nv(G))
    for i in 1:length(P)
        C[P[i]] .= C_q[i] ./ length(P[i])
    end
    return C, P
end

function quotient_graph(G::AbstractGraph{T};
    weights::Union{SparseMatrixCSC{<:Number,Int},Nothing}=nothing, args...) where {T}
    P = refine_fixpoint(G; weights=weights, args...)

    k = length(P)
    G₀ = SimpleDiGraph(k)

    𝜑 = Dict{T,Color}()
    sizehint!(𝜑, nv(G))
    for (color, nodes) in enumerate(P)
        for x in nodes
            𝜑[x] = color
        end
    end

    if weights === nothing
        weights::SparseMatrixCSC{Float64,Int} = adjacency_matrix(G)
    end

    P_sparse = partition_matrix(P)
    neighbor = weights * P_sparse

    @assert sum(neighbor, dims=2) == sum(weights, dims=2)

    C₀::Matrix{Float64} = zeros(k, k)
    for i in eachindex(P)
        X::Vector{Int64} = P[i]
        C₀[i, :] = sum(neighbor[X, :], dims=1)
    end

    # @assert sum(C₀) <= sum(weights)

    for i in findall(w -> w != 0.0, C₀)
        u, v = i[1], i[2]
        add_edge!(G₀, u, v)
    end

    return G₀, C₀, 𝜑, P
end

function quotient_centrality_silly(G, colors)
    P = refine_fixpoint(G; early_stop=colors)
    C₀::Vector{Float64} = zeros(nv(G))
    for Pᵢ in P
        C₀ .+= betweenness_centrality(G, [Pᵢ[1]], normalize=false) * length(Pᵢ)
    end
    return C₀
end

function quotient_centrality(G_q, W, P, G)
    n = nv(G_q)
    C = zeros(n)
    all_paths = 0

    L = [length(X) for X in P]

    for u in 1:n
        P_deg = zeros(n)
        # there is one path from u to u, the empty set
        P_deg[u] = 1

        ds = dijkstra_shortest_paths(G_q, u, trackvertices=true, allpaths=true)

        for v in ds.closest_vertices
            if ds.pathcounts[v] == 0
                break
            end

            if v != u
                pred = ds.predecessors[v]
                paths = sum(P_deg[pred] .* W[pred, v] ./ max.(L[v], L[pred]))
                @assert paths > 0
                P_deg[v] = paths
                all_paths += paths
            end
        end

        P_count = zeros(n)

        for v in reverse(ds.closest_vertices)
            if ds.pathcounts[v] == 0 || ds.dists[v] <= 1
                continue
            end

            if v <= u
                continue
            end

            pred = ds.predecessors[v]
            if pred != [u]
                P_count[pred] .+= (P_deg[v] + L[v] * P_count[v]) .* (P_deg[pred] ./ sum(P_deg[pred]))
            end
        end

        C .+= L[u] * P_count

        # if L[u] > 1 && W[u, u] == 0.0
        #     # expand super-node self-loops
        #     w = neighbors(G_q, u)
        #     d = minimum(ds.dists[w])
        #     close = filter(x -> ds.dists[x] == d, w)
        #     degs = W[close, u] ./ L[close]
        #     if any(degs .> 1.0)
        #         C[close] .+= binomial(L[u], 2) * (degs .- 1) ./ sum(degs .- 1)
        #     end
        # end
        if L[u] > 1 && W[u, u] == 0.0
            C_loop = betweenness_centrality(G, [P[u][1], P[u][2]]; normalize=false)
            C .+= binomial(L[u], 2) * [mean(C_loop[P[i]]) for i in 1:n]
        end
    end

    return C
end
