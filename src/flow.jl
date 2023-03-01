"""Functions related to computing flows on network graphs."""
module Flow

export lifted_maxflow

using QuasiStableColors: q_color, super_graph, node_map
using Graphs, SparseArrays

"""
    lifted_maxflow(
        solver,
        G::Graph,
        s::Int,
        t::Int,
        q = 0.0,
        n_colors = Inf,
        weights::SparseMatrixCSC{<:Number,Int} = nothing,
    )

Return an equivalent approximate maximum flow problem from `s` to `t` in flow network `G`. Capacities
defined using `weights`; if none provided, unit capacities assumed.
Uses a quasi-stable coloring with maximum error `q` or `n_colors` colors, whichever is
smaller. `Solver` should be `GraphsFlows.maximum_flow` or similiar function.

Example of computing the approximate maximum flow:
    using GraphsFlows: maximum_flow
    f_approx, _ = lifted_maxflow(maximum_flow, G, s, t, C, n_colors=25)
"""
function lifted_maxflow(solver, G, s::Int, t::Int; args...)::Number
    C = q_color(G; special=Set([s, t]), args...)
    G₀, C₀ = super_graph(C)
    𝜑 = node_map(C)
    s₀, t₀ = 𝜑[s], 𝜑[t]
    f₀, _ = solver(G₀, s₀, t₀, C₀)
    return f₀
end
end
