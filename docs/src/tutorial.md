# Tutorial

Let's explore how to use the QuasiStableColors library with an example.


## Installation
First, install the library. One option is to run the following code:
```@example
using Pkg
Pkg.add("QuasiStableColors")
```
*alternatively*, open the Julia REPL and press `]` to open the package manager.
Then type `install QuasiStableColors`:
```julia
 pkg> install QuasiStableColors
```

After installing either way, now activate the libary:
```@example
using QuasiStableColors
```
## Coloring Graphs
Let's load an example graph.

```@example
using Graphs
n = 8
g = Graphs.SimpleGraphs.dorogovtsev_mendes(n, seed=0)

using GraphMakie, CairoMakie
f, ax, p = graphplot(g)
hidedecorations!(ax); hidespines!(ax) #hide
save("graph.svg", f); nothing #hide
```

```@raw html
<img src="/graph.svg" alt="Illustration of network graph">
```
