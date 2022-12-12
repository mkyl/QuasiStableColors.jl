Experiment Scripts
======

From the paper "Quasi-stable Coloring for Graph Compression: Approximating Max-Flow,
Linear Optimization and Centrality" by Moe Kayali and Dan Suciu, appearing in VLDB '23.

Contents:

- `experiments/`
  - `figures.jl` color to make the karate club figure 
  - `histograms.jl` builds the histograms of color sizes
  - `measure-centrality.jl` accuracy-speedup trade-off curves for centrality
  - `measure-flow.jl` accuracy-speedup trade-off curves for max-flow
  - `measure-lp.jl` accuracy-speedup trade-off curves for linear programs
  - `plot.jl` draws accuracy-speedup curves from the three previous experiments
  - `robustness.jl` resilience to random pertubation of edges experiment 
  - `scale.jl` measuring the number of color for a particular q-error and ε-error
  - `runtime.jl` draw plots that break down runtime into various steps
