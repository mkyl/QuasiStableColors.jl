"""Result of quasi-stable coloring algorithm."""
struct Coloring{T<:Integer}
    partitions::Vector{Vector{T}}
    max_q::Number
    stats_out
    stats_in
end

partitions(C::Coloring) = C.partitions

max_q_err(C::Coloring) = C.max_q

"""Return a `Dict` mapping from node id to its color."""
function node_map(C::Coloring{T}) where {T<:Integer}
    𝜑 = Dict{T,Color}()
    # sizehint!(𝜑, nv(G))
    P = partitions(C)
    for (color, nodes) in enumerate(P)
        for x in nodes
            𝜑[x] = color
        end
    end
    return 𝜑
end


"""Build a lifted graph where each node is color. Also called a quotient graph."""
function super_graph(C::Coloring{T}; agg=sum) where {T<:Integer}
    P = partitions(C)
    k = length(C.partitions)

    weights = zeros(k, k)
    for i in eachindex(P)
        X::Vector{Int64} = P[i]
        weights[i, :] .= transpose(sum(C.stats_out.neighbor[X, :], dims=1))
    end

    G = SimpleDiGraph(weights)

    return G, weights
end
