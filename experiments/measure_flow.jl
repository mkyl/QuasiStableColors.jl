using Graphs
using GraphsFlows

using CSV
using DataFrames

include("../flow.jl")
include("../datasets.jl")

ENV["JULIA_DEBUG"] = Main

DATASETS = ["sawtooth0", "cells", "sawtooth1", "tsukuba0", "tsukuba1", "tsukuba2", "venus0", "venus1", "venus2", "simcells",]

function measure_flow()
    parent = pwd()
    mkpath(joinpath(parent, "results"))
    results = joinpath(parent, "Lifted Inference", "results", "flow.tsv")
    if isfile(results)
        df = CSV.read(results, DataFrame; delim='\t')
    else
        df = DataFrame(name=String[], colors=Int[], time=Float64[], size=Int[],
            objective=Float64[])
    end

    compile = false
    for x in DATASETS
        G, C, s, t = getfield(datasets, Symbol(x))()

        # once for JIT compiler
        if !compile
            lifted_maxflow(G, s, t; weights=C, early_stop=4)
            compile = true
        end

        if !any((df.name .== x) .& (df.colors .== 0))
            stats = @timed maximum_flow(G, s, t, C)
            obj, _ = stats.value
            time = stats.time
            size = nv(G)
            push!(df, [x 0 time size obj])
            CSV.write(results, df; delim='\t')
        end

        for i in 3:2:35
            if !any((df.name .== x) .& (df.colors .== i))
                stats = @timed lifted_maxflow(G, s, t; weights=C, early_stop=i)
                obj = stats.value
                time = stats.time
                push!(df, [x i time i obj])
                CSV.write(results, df; delim='\t')
            end
        end
    end

    CSV.write(results, df; delim='\t')
end

# measure_flow()
