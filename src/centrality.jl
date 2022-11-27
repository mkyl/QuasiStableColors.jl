"""Functions for approximating node centrality."""
module Centrality

export approx_betweenness_centrality

using QuasiStableColors
using Graphs

"""Approximate betweenness centrality using a size `colors` q-stable coloring."""
function approx_betweenness_centrality(G; args...)
    P = q_color(G; args...)
    C₀::Vector{Float64} = zeros(nv(G))
    for Pᵢ in P
        C₀ .+= betweenness_centrality(G, [Pᵢ[1]], normalize=false) * length(Pᵢ)
    end
    C₀ = C₀ ./ nv(G)
    return C₀
end
end
