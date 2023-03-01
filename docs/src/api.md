# API
This is a comprehensive list functions implemented in this library. For a guided introduction to these methods, check out [the tutorial](@ref Tutorial) first.
### Page contents
```@contents
Pages = ["api.md"]
```
## Coloring
These are the core methods to compute quasi-stable colorings of graphs.
```@docs
q_color(G;
    weights=nothing,
    special::Set{T}=Set{T}(),
    warm_start::Vector{Vector{T}}=Vector{Vector{T}}(),
    early_stop=Inf,
    eps::Float64=0.0) where {T}
```
Coloring an example graph:
```@example
using QuasiStableColors: q_color
using Graphs

# Make a chain graph
G = SimpleGraph(4)
add_edge!(G, 1, 1)
add_edge!(G, 1, 2)
add_edge!(G, 2, 3)
add_edge!(G, 3, 4)
add_edge!(G, 4, 4)

# Color it
C = q_color(G, q=2.0)
```

```@docs
Coloring
node_map
super_graph
```

For the specific case of bipartite graphs, consider:
```@docs
refine_bipartite
```

Finally, for comparison with prior work, we provide the following method.
```@docs
refine_stable(G)
```
## Maximum-flow
This method for approximating [maximum-flow](https://en.wikipedia.org/wiki/Maximum_flow_problem) in a flow network.
```@docs
QuasiStableColors.Flow.lifted_maxflow
```
An example of how to compute an approximate maximum flow:
```@example
using QuasiStableColors.Flow: lifted_maxflow
using Graphs
using GraphsFlows: maximum_flow

E = [Edge(1,2), Edge(1, 3), Edge(1, 4), Edge(2, 5), Edge(3, 5), Edge(4, 5)]
G = SimpleGraphFromIterator(E)

lifted_maxflow(maximum_flow, G, 1, 5; q=1.0)
```

## Betweenness centrality
This method is provided for approximating [betweenness centrality](https://en.wikipedia.org/wiki/Betweenness_centrality).
```@docs
QuasiStableColors.Centrality.approx_betweenness_centrality
```
To compute the estimated betweenness centralities on a sample graph:
```@example
using Graphs
using QuasiStableColors.Centrality: approx_betweenness_centrality

# Construct example graph
E = [Edge(1, 5), Edge(2, 4), Edge(2, 5), Edge(2, 8), Edge(3, 5), Edge(3, 9),
    Edge(6, 9), Edge(7, 8), Edge(8, 9)]
G = SimpleGraphFromIterator(E)

# Compute approximate centrality
C₀ = approx_betweenness_centrality(G, q=0.0)
println("Approximate centrality: $C₀")

# Compute exact centrality for comparision
using Graphs: betweenness_centrality
C = betweenness_centrality(G, normalize=false)
println("      Exact centrality: $C")
```
## Linear programming
This section details methods related to approximating linear optimization.
```@docs
QuasiStableColors.Optimize.lifted_minimize
QuasiStableColors.Optimize.lifted_maximize
```
Consider the following example linear system from [our paper](https://arxiv.org/abs/2211.11912):
```math
\begin{align*}
      \text{maximize } & 9x_1+10x_2+50x_3 \\
      \text{where } & 4x_1+8x_2+2x_3 \leq 20 \\
                       & 6x_1 + 5x_2 + x_3 \leq 20 \\
                       & 7x_1 + 4x_2 + 2x_3 \leq 21 \\
                       & 3x_1 + x_2 + 22x_3 \leq 50 \\
                       & 2x_1 + 3x_2 + 21x_3 \leq 51
\end{align*}
```
This system has optimum value ``c^T x^\ast = 128.2``. Let's compute the approximate minimum:
```@example
using QuasiStableColors.Optimize: lifted_maximize
using Tulip

solver = Tulip.Optimizer()

A = [4.0 8 2; 6 5 1; 7 4 2; 3 1 22; 2 3 21]
b = [20.0, 20, 21, 50, 51]
c = [9.0, 10, 50]
z = lifted_maximize(solver, A, b, c; q=5.0)
"estimated value: $z"
```
This gives us an estimated value of ``c^T x^\ast`` within ``1\%`` of the true value.
