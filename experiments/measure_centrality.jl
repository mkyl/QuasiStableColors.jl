using QPSReader
using SparseArrays
using CSV
using DataFrames

include("../centrality.jl")
include("../datasets.jl")

DATASETS = [datasets.astrophysics, datasets.enron, datasets.facebook,
    datasets.epinions, datasets.deezer]

ENV["JULIA_DEBUG"] = Main

function measure_centrality()
    parent = pwd()
    mkpath(joinpath(parent, "results"))
    results = joinpath(parent, "Lifted Inference", "results", "centrality.tsv")
    if isfile(results)
        df = CSV.read(results, DataFrame; delim='\t')
    else
        df = DataFrame(name=String[], colors=Int[], time=Float64[],
            corr=Float64[])
    end

    for x in DATASETS
        G = x()

        stats = @timed betweenness_centrality(G, normalize=false)
        C = stats.value
        time = stats.time

        # for the JIT
        quotient_centrality_silly(G, 5)
        for i in [5, 25, 50, 100, 150, 200]
            if !any((df.name .== x) .& (df.colors .== i))
                stats = @timed quotient_centrality_silly(G, i)
                C₀ = stats.value
                time₀ = stats.time
                runtime = time₀ / time * 100.0
                corr = corspearman(C, C₀)
                name = String(Symbol(x))
                push!(df, [name i runtime corr])
            end
        end
        CSV.write(results, df; delim='\t')
    end
end

measure_centrality()
