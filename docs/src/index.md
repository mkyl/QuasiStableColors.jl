# Quasi-Stable Coloring
*Graph compression for performant approximations*

QuasiStableColors.jl is a library that enables compressing graphs. The compression is
lossy, meaning that not all edges are preserved. The compressed graphs are computed
using an algorithm called *quasi-stable coloring*, introduced in the paper...

## Citation Format
If you use this algorithm, we ask that you cite our paper:
```
@misc {
    title = "Quasi-stable Coloring for Graph Compression: Approximating Max-Flow, Linear Programs, and Centrality"
    author = "Moe Kayali and Dan Suciu"
}
```

## Applications
Developed applications for quasi-stable coloring include the maximum-flow/minimum-cut
problems, betweenness centrality computation and linear optimization.  

## Getting Started
The [tutorial section](tutorial) contains all the code you need to get started with
using this library.

## API
Reference documentation for all the public functions of this library can be found in the
[API](api/) section.

## Researchers' Guide
Want to extend quasi-stable coloring to a new domain? Perhaps you want to develop a
variant of approximate colorings? The [internals](internals) section
covers... 
