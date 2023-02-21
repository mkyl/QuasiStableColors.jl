"""Result of quasi-stable coloring algorithm."""
struct qColoring{T<:Integer}
    partitions::Vector{Vector{T}}
    max_q::Number
    k::Integer
end

partitions(C::qColoring) = C.partitions

max_q_err(C::qColoring) = C.max_q

"""Return a `Dict` mapping from node id to its color."""
function node_map(C::qColoring)
    𝜑 = Dict{T,Color}()
    sizehint!(𝜑, nv(G))
    for (color, nodes) in enumerate(P)
        for x in nodes
            𝜑[x] = color
        end
    end
    return 𝜑
end


"""Build a lifted graph where each node is color. Also called a quotient graph."""
function super_graph(C::qColoring; agg=sum)
    k = length(C.partitions)
    G = SimpleDiGraph(k)
    return G
end
