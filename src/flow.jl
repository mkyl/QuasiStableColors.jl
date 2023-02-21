"""Functions related to computing flows on network graphs."""
module Flow

export lifted_maxflow

using QuasiStableColors
using Graphs, GraphsFlows, SparseArrays

include("api.jl")

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
    C = q_color(G; adjacency=weights, special=Set([s, t]), args...)
    G₀, C₀ = super_graph(C)
    𝜑 = node_map(C)
    s₀, t₀ = 𝜑[s], 𝜑[t]
    f₀, _ = maximum_flow(G₀, s₀, t₀, C₀)
    return f₀
end
end
