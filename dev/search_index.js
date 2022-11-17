var documenterSearchIndex = {"docs":
[{"location":"api/#API","page":"API","title":"API","text":"","category":"section"},{"location":"api/#Coloring","page":"API","title":"Coloring","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"refine_fixpoint(G;\n    weights=nothing,\n    special::Set{T}=Set{T}(),\n    warm_start::Vector{Vector{T}}=Vector{Vector{T}}(),\n    early_stop=Inf,\n    eps::Float64=0.0) where {T}\nrefine_stable(G)\nrefine_bipartite","category":"page"},{"location":"api/#QuasiStableColors.refine_fixpoint-Union{Tuple{Any}, Tuple{T}} where T","page":"API","title":"QuasiStableColors.refine_fixpoint","text":"refine_fixpoint(\n    G::AbstractGraph{T},\n    weights = nothing,\n    special =  Set{T}(),\n    warm_start = Vector{Vector{T}}(),\n    early_stop = Inf,\n    eps = 0.0,\n)\n\nCompute a quasi-stable coloring for the undirected graph G. Typically, you should  set one of:\n\neps: maximum q-error allowed\nearly_stop: number of colors to use\n\nOptional parameters:\n\nwarm_start: coloring to refine. If not provided, start using trivial (single color) partitioning assumed.\n\n\n\n\n\n","category":"method"},{"location":"api/#QuasiStableColors.refine_stable-Tuple{Any}","page":"API","title":"QuasiStableColors.refine_stable","text":"refine_stable(G::AbstractGraph{T})\n\nCompute the stable coloring for the undirected graph G. Provided for comparasion.\n\n\n\n\n\n","category":"method"},{"location":"api/#QuasiStableColors.refine_bipartite","page":"API","title":"QuasiStableColors.refine_bipartite","text":"Equivalent to refine_fixpoint but optimized for bipartite graphs. Faster but less general.\n\n\n\n\n\n","category":"function"},{"location":"api/#Maximum-flow","page":"API","title":"Maximum-flow","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"QuasiStableColors.Flow.lifted_maxflow","category":"page"},{"location":"api/#QuasiStableColors.Flow.lifted_maxflow","page":"API","title":"QuasiStableColors.Flow.lifted_maxflow","text":"Compute the maximum flow from s to t in flow network G. Capacities defined using weights=.\n\n\n\n\n\n","category":"function"},{"location":"api/#Betweenness-centrality","page":"API","title":"Betweenness centrality","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"QuasiStableColors.Centrality.approx_betweenness_centrality","category":"page"},{"location":"api/#QuasiStableColors.Centrality.approx_betweenness_centrality","page":"API","title":"QuasiStableColors.Centrality.approx_betweenness_centrality","text":"Approximate betweenness centrality using a size colors q-stable coloring.\n\n\n\n\n\n","category":"function"},{"location":"api/#Linear-programming","page":"API","title":"Linear programming","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"QuasiStableColors.Optimize.lifted_minimize","category":"page"},{"location":"api/#QuasiStableColors.Optimize.lifted_minimize","page":"API","title":"QuasiStableColors.Optimize.lifted_minimize","text":"Approximate linear program min c^T x where A x >= b, x >=0.\n\n\n\n\n\n","category":"function"},{"location":"applications/#Applications","page":"Applications","title":"Applications","text":"","category":"section"},{"location":"applications/#Max-flow/min-cut","page":"Applications","title":"Max-flow/min-cut","text":"","category":"section"},{"location":"applications/#Betweenness-Centrality","page":"Applications","title":"Betweenness Centrality","text":"","category":"section"},{"location":"applications/#Linear-Optimization","page":"Applications","title":"Linear Optimization","text":"","category":"section"},{"location":"#Quasi-Stable-Coloring","page":"Quasi-Stable Coloring","title":"Quasi-Stable Coloring","text":"","category":"section"},{"location":"","page":"Quasi-Stable Coloring","title":"Quasi-Stable Coloring","text":"Graph compression for performant approximations","category":"page"},{"location":"","page":"Quasi-Stable Coloring","title":"Quasi-Stable Coloring","text":"QuasiStableColors.jl is a library for compressing graphs and approximating graph algorithms. The compressed graphs are computed using an algorithm called quasi-stable coloring, which results in a much smaller graph while preserving its key properties. This approach is introduced in the research paper \"Quasi-stable Coloring for Graph Compression.\"","category":"page"},{"location":"","page":"Quasi-Stable Coloring","title":"Quasi-Stable Coloring","text":"A major advantage of this approach is that many algorithms can be computed directly on the compressed graph, without needing decompression. This results in an effective  approximation of many graph algorithms. Applications implemented in this library are:","category":"page"},{"location":"","page":"Quasi-Stable Coloring","title":"Quasi-Stable Coloring","text":"Betweenness centrality\nMaximum-flow/minimum-cut\nLinear optimization","category":"page"},{"location":"#Citation-Format","page":"Quasi-Stable Coloring","title":"Citation Format","text":"","category":"section"},{"location":"","page":"Quasi-Stable Coloring","title":"Quasi-Stable Coloring","text":"If you use this library, we ask that you cite our paper:","category":"page"},{"location":"","page":"Quasi-Stable Coloring","title":"Quasi-Stable Coloring","text":"@article {\n    title = {Quasi-stable Coloring for Graph Compression: Approximating Max-Flow, Linear Programs, and Centrality},\n    author = {Moe Kayali and Dan Suciu},\n    journal = {Proc. VLDB Endow.},\n    volume = {15},\n    year = {2022},\n}","category":"page"},{"location":"","page":"Quasi-Stable Coloring","title":"Quasi-Stable Coloring","text":"Quasi-stabling coloring was developed at the School of Computer Science in the University of Washington, Seattle.","category":"page"},{"location":"#Resources","page":"Quasi-Stable Coloring","title":"Resources","text":"","category":"section"},{"location":"","page":"Quasi-Stable Coloring","title":"Quasi-Stable Coloring","text":"Tutorial: examples and code to get started with using this library.\nQuestions: ask questions using the quasi-stable-colors tag on StackOverflow.\nReference: Reference documentation for all the public functions of this library can be found in the API section. Specific sections explain the applications: maximum-flow/minimum-cut problems, betweenness centrality computation and linear optimization.\nResearchers' Guide: Want to extend quasi-stable coloring to a new domain? Perhaps you want to develop a variant of approximate colorings? The Internals section covers these topics.","category":"page"},{"location":"internals/#Internals","page":"Internals","title":"Internals","text":"","category":"section"},{"location":"internals/","page":"Internals","title":"Internals","text":"This package implements the Rothko algorithm, presented in the paper.","category":"page"},{"location":"tutorial/#Tutorial","page":"Tutorial","title":"Tutorial","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Hello","category":"page"}]
}
