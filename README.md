# QuasiStableColors.jl
*Graph compression for performant approximations*
![Project logo](docs/src/assets/logo.png)

| Documentation  | Build Status |
| ------------- | ------------- |
| [Stable](https://mkyl.github.io/QuasiStableColors.jl/stable/)  | ![continuous integration status badge](https://github.com/mkyl/QuasiStableColors.jl/actions/workflows/CI.yml/badge.svg)  |


Code and experiments for the upcoming paper *"Quasi-stable Coloring
for Graph Compression: Approximating Max-Flow, Linear Optimization
and Centrality"* by Moe Kayali and Dan Suciu, to appear in VLDB '23. Find the paper [here](https://arxiv.org/abs/2211.11912).

Contents:
- 'src/'
  - `QuasiStableColors.jl`: implementation of the Rothko algorithm
  - `datasets.jl`: handlers to ingest datasets
  - `flow.jl`: code for max-flow approximation
  - `optimize.jl`: code for linear program approximation
  - `centrality.jl`: code for betweenness centrality approximation
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
-  `tests/` correctness checks for the implementation code
