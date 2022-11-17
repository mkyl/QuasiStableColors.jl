# API

## Coloring
```@docs
refine_fixpoint(G;
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
QuasiStableColors.Flow.lifted_maxflow
```

## Betweenness centrality
```@docs
QuasiStableColors.Centrality.approx_betweenness_centrality
```


## Linear programming
```@docs
QuasiStableColors.Optimize.lifted_minimize
```
