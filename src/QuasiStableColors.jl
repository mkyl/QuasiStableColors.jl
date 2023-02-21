"""Compute quasi-stable colorings for graphs. Includes approximating max-flow, linear
programs and betweenness centrality."""
module QuasiStableColors

include("api.jl")
export qColoring, Color, partition, node_map, max_q_err, super_graph

include("refine.jl")
export refine_stable, q_color

include("bipartite.jl")
export refine_bipartite

include("centrality.jl")
export approx_betweenness_centrality

include("flow.jl")
export lifted_maxflow

include("optimize.jl")
export lifted_maximize, lifted_minimize, maximize, minimize

end
