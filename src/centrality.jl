"""Functions for approximating node centrality."""
module Centrality

export approx_betweenness_centrality

using QuasiStableColors: q_color, partitions
using Graphs

"""
    approx_betweenness_centrality(
        G::Graph,
        q::Number,
        n_colors::Int,
    )

Approximate betweenness centrality using a q-stable coloring with maximum error `q` or
size `n_colors`, whichever is smaller."""
function approx_betweenness_centrality(G; args...)
    colors = q_color(G; args...)
    P = partitions(colors)
    C₀::Vector{Float64} = zeros(nv(G))
    for Pᵢ in P
        C₀ .+= betweenness_centrality(G, [Pᵢ[1]], normalize=false) * length(Pᵢ)
    end
    C₀ = C₀ ./ nv(G)
    return C₀
end
end
