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
refine_bipartite
QuasiStableColoring
```

Finally, for comparison with prior work, we provide:
```@docs
refine_stable(G)
```
## Maximum-flow
```@docs
QuasiStableColors.Flow
QuasiStableColors.Flow.lifted_maxflow
```
An example of how to compute an approximate maximum flow:

## Betweenness centrality
This method is provided for approximating [betweenness centrality](https://en.wikipedia.org/wiki/Betweenness_centrality).
```@docs
QuasiStableColors.Centrality.approx_betweenness_centrality
```
To compute the estimated betweenness centralities on a sample graph:
```@example
1+1
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

A = [4 8 2; 6 5 1; 7 4 2; 3 1 22; 2 3 21]
b = [20, 20, 21, 50, 51]
c = [9, 10, 50]
z = lifted_maximize(A, b, c; q=5.0)
"estimated value: $z"
```
This gives us an estimated value of ``c^T x^\ast`` within ``1\%`` of the true value.
