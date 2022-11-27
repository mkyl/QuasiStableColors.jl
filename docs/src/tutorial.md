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
After installing activate the libary:
```@example coloring
using QuasiStableColors
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

Now, let's generate a quasi-stable coloring `C` where we allow at most one edge
error (*i.e.* $q=1$).
```@example coloring
C = QuasiStableColors.refine_fixpoint(g, eps=1.0)
```

We get a four-partition coloring. Let's assign a unique graphical color to each:

```@example coloring
using Colors
palette = distinguishable_colors(length(C))
palette = distinguishable_colors(length(C), [RGB(1,1,1), RGB(0,0,0)], dropseed=true) #hide
color_map = Array{Colorant}(undef, nv(g))
for (i, x) in enumerate(C)
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

## Approximation 
Now that we have an idea of how coloring works, we examine why it's 
useful as an approximation. We are going to try to compute maximum-flow on the
graph.

First, let's define the [maximum flow problem]
(https://en.wikipedia.org/wiki/Maximum_flow_problem). Pick two special nodes and
call them *source* and *sink*. Now imagine each edge is a pipe of water. The
question max-flow problem asks is: *what is the largest volume of water can we
simultaneously push from the source to the sink?*  For simplicity we assume each
pipe has the same size (aka *unit capacities*), carrying at most a gallon per
minute.

For example, consider the flows between nodes $X$ and $Y$. One valid way to push
water, called a *flow*, is to send one gallon / min of water down to Z and another
gallon down to W. A total of two gallons would arrive at $Y$. But is the most water
we could send?

Many algorithms have been devised to solve this classic problem. Let's consider how
we can use the lifted graph.

If we run the algorithm on the ground graph, we find in fact this is the correct
solution:


