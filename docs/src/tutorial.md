# Tutorial

Let's explore how to use the QuasiStableColors library with an example.

```@setup libraries
using QuasiStableColors
using Graphs
```


## Installation
First, install the library. One option is to run the following code:
```julia
using Pkg
Pkg.add("QuasiStableColors")
```
After installing activate the library:
```@example coloring
using QuasiStableColors
```

Let's use `QSC` as shorthand for the library name:
```@example coloring
const QSC = QuasiStableColors
nothing #hide
```
## Coloring Graphs
Let's create a simple example graph.

```@example coloring
using Graphs
n = 8
g = Graphs.SimpleGraphs.dorogovtsev_mendes(n, seed=0)

using GraphMakie, CairoMakie
graphplot(g)
f, ax, p = graphplot(g, node_size=32) #hide
hidedecorations!(ax); hidespines!(ax) #hide
save("graph.svg", f); nothing #hide
```

![Example network graph from above code](graph.svg)

Now, let's generate a quasi-stable coloring `C` where we allow at most one edge error (*i.e.* $q=1$).
```@example coloring
C = QSC.q_color(g, q=1.0)
```

We get a four-partition coloring. How are the nodes grouped?
```@example coloring
P = QSC.partitions(C)
```
Let's assign a unique visual color to each group:
```@example coloring
using Colors
palette = distinguishable_colors(length(P))
palette = distinguishable_colors(length(P), [RGB(1,1,1), RGB(0,0,0)], dropseed=true) #hide
color_map = Array{Colorant}(undef, nv(g))
for (i, x) in enumerate(P)
    color_map[x] .= palette[i]
end 
```
and draw the graph again:
```@example coloring
graphplot(g, node_color=color_map)
f_c, ax_c, _ = graphplot(g, node_color=color_map, node_size=32) #hide
hidedecorations!(ax_c); hidespines!(ax_c) #hide
save("graph-colors.svg", f_c); nothing #hide
```

![Same graph, now colored according to above partition](graph-colors.svg)

## Summary Graphs
Let's build a summary graph using the coloring `C`. Each supernode represents one color:

```@example coloring
G_super, weights = QSC.super_graph(C)
f, ax, p = graphplot(G_super, node_color=palette, node_size=32) #hide
hidedecorations!(ax); hidespines!(ax) #hide
save("super-graph.svg", f); nothing #hide
```
![Graph of supernodes implied by coloring](super-graph.svg)

## Approximation
This section on how to use the coloring for graph algorithms is upcoming--for now see
the [maximum-flow section](@ref Maximum-flow) of the API reference as an example.
