"""Compute quasi-stable colorings for graphs. Includes approximating max-flow, linear
programs and betweenness centrality."""
module QuasiStableColors

include("refine.jl")
include("bipartite.jl")
include("centrality.jl")
include("flow.jl")
include("optimize.jl")

export refine_stable, q_color, refine_bipartite, Color, q_color_dir
export approx_betweenness_centrality
export lifted_maxflow
export lifted_maximize, lifted_minimize, maximize, minimize
export QuasiStableColoring
end
