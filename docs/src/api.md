# API

## Coloring
```@docs
q_color(G;
    weights=nothing,
    special::Set{T}=Set{T}(),
    warm_start::Vector{Vector{T}}=Vector{Vector{T}}(),
    early_stop=Inf,
    eps::Float64=0.0) where {T}
refine_stable(G)
refine_bipartite
```

## Maximum-flow
```@docs
QuasiStableColors.Flow
QuasiStableColors.Flow.lifted_maxflow
```

## Betweenness centrality
```@docs
QuasiStableColors.Centrality
QuasiStableColors.Centrality.approx_betweenness_centrality
```


## Linear programming
```@docs
QuasiStableColors.Optimize
QuasiStableColors.Optimize.lifted_minimize
QuasiStableColors.Optimize.lifted_maximize
```
