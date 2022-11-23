# Quasi-Stable Coloring
*Graph compression for performant approximations*

QuasiStableColors.jl is a library for compressing graphs and approximating graph
algorithms. The compressed graphs are computed using an algorithm called *quasi-stable
coloring*, which results in a much smaller graph while preserving its key properties.
This approach is introduced in the research paper "Quasi-stable Coloring for Graph
Compression," which can be read [here](https://arxiv.org/abs/2211.11912).

A major advantage of this approach is that many algorithms can be computed directly on
the compressed graph, *without needing decompression*. This results in an effective 
approximation of many graph algorithms. Applications implemented in this library are:
 - Betweenness centrality
 - Maximum-flow/minimum-cut
 - Linear optimization
## Citation Format
If you use this library, we ask that you cite our paper:
```
@article {
    title = {Quasi-stable Coloring for Graph Compression: Approximating Max-Flow, Linear Programs, and Centrality},
    author = {Moe Kayali and Dan Suciu},
    journal = {Proc. VLDB Endow.},
    volume = {15},
    year = {2022},
}
```
Quasi-stabling coloring was developed in the [School of Computer
Science](https://www.cs.washington.edu) at the University of Washington, Seattle.
## Resources
- **[Tutorial](@ref)**: examples and code to get started with using this library.
- **[Reference](@ref API)**: Reference documentation for all the public functions of this library can be found in the [API](@ref) section. Specific sections explain the applications: maximum-flow/minimum-cut problems, betweenness centrality computation and linear optimization.
- **Questions**: ask questions and report bugs using [Github issues](https://github.com/mkyl/QuasiStableColors.jl/issues).
- **[Researchers' Guide](@ref Internals)**: Want to extend quasi-stable coloring to a new domain? Perhaps you want to develop a variant of approximate colorings? The [Internals](@ref) section covers these topics.
