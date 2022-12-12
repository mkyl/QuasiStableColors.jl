"""Functions related to computing flows on network graphs."""
module Flow

export lifted_maxflow

using QuasiStableColors
using Graphs, GraphsFlows, SparseArrays

include("misc.jl")

function _quotient_graph(G::AbstractGraph{T};
    weights::Union{SparseMatrixCSC{<:Number,Int},Nothing}=nothing, args...) where {T}
    P = q_color(G; weights=weights, args...)

    k = length(P)
    G₀ = SimpleDiGraph(k)

    𝜑 = Dict{T,Color}()
    sizehint!(𝜑, nv(G))
    for (color, nodes) in enumerate(P)
        for x in nodes
            𝜑[x] = color
        end
    end

    P_sparse = partition_matrix(P)

    if weights === nothing
        @debug "Assuming unit capacities"
        weights::SparseMatrixCSC{Float64,Int} = adjacency_matrix(G)
        weights.nzval .= 1.0
    end

    neighbor = weights * P_sparse

    @assert sum(neighbor, dims=2) == sum(weights, dims=2)

    C₀::Matrix{Float64} = zeros(k, k)
    for i in eachindex(P)
        X::Vector{Int64} = P[i]
        C₀[i, :] = sum(neighbor[X, :], dims=1)
    end

    for i in findall(w -> w != 0.0, C₀)
        u, v = i[1], i[2]
        add_edge!(G₀, u, v)
    end

    return G₀, C₀, 𝜑, P
end

"""
    lifted_maxflow(
        G::Graph,
        s::Int,
        t::Int,
        q = 0.0,
        n_colors = Inf,
        weights::SparseMatrixCSC{<:Number,Int} = nothing,
    )

Compute the approximate maximum flow from `s` to `t` in flow network `G`. Capacities
defined using `weights`; if none provided, unit capacities assumed.
Uses a quasi-stable coloring with maximum error `q` or `n_colors` colors, whichever is
smaller."""
function lifted_maxflow(G, s::Int, t::Int; args...)::Number
    G₀, C₀, 𝜑, _ = _quotient_graph(G; special=Set([s, t]), args...)
    s₀, t₀ = 𝜑[s], 𝜑[t]
    f₀, _ = maximum_flow(G₀, s₀, t₀, C₀)
    return f₀
end
end
