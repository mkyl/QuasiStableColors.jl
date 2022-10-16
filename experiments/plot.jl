using CSV
using Plots
using DataFrames
using Statistics: mean

DATASETS = ["qap15.mps", "ex10.mps", "supportcase10.mps", "nug08-3rd.mps",]
DATASETS_FLOW = ["simcells", "cells", "tsukuba0", "tsukuba2", "venus0", "venus1", "sawtooth0", "sawtooth1"]

movingaverage(g, n) = [i < n ? mean(g[begin:i]) : mean(g[i-n+1:i]) for i in 1:length(g)]

function make_lp_plot()
    plot()
    results = joinpath(pwd(), "Lifted Inference", "results", "lp.tsv")
    df = CSV.read(results, DataFrame; delim='\t')

    for x in DATASETS
        df_x = subset(df, :name => ByRow(==(x)))
        baseline = subset(df_x, :colors => ByRow(==(0)))
        baseline = baseline[1, :]
        df_x = subset(df_x, :colors => ByRow(!=(0)))

        sort!(df_x, [:time])
        df_x.time = df_x.time / baseline.time * 100
        df_x.objective .= (df_x.objective ./ baseline.objective)
        df_x.objective .= max.(df_x.objective, 1 ./ df_x.objective)
        plot!(df_x.time, df_x.objective, label=splitext(x)[1], lw=2, markershape=:auto)
    end

    hline!([1.0], label="exact", linestyle=:dash, lw=2)
    ylabel!("relative error")
    xlabel!("Runtime (% of baseline)")

    plot!(size=(400, 300), ylim=(0.95, 2.35), xlim=(0.0, 1.0), thickness_scaling=1.25,
        legend=:topright, left_margin=-3Plots.mm, bottom_margin=-3Plots.mm)
    savefig("lp.pdf")
end

function make_lp_color_plot()
    plot()
    results = joinpath(pwd(), "Lifted Inference", "results", "lp.tsv")
    df = CSV.read(results, DataFrame; delim='\t')

    for x in DATASETS
        df_x = subset(df, :name => ByRow(==(x)))
        baseline = subset(df_x, :colors => ByRow(==(0)))
        baseline = baseline[1, :]
        df_x = subset(df_x, :colors => ByRow(!=(0)))
        sort!(df_x, [:colors])
        df_x.time = df_x.time / baseline.time * 100
        df_x.objective .= (df_x.objective ./ baseline.objective)
        df_x.objective .= max.(df_x.objective, 1 ./ df_x.objective)
        plot!(df_x.colors, df_x.objective, label=splitext(x)[1], lw=2, markershape=:auto)
    end

    hline!([1.0], label="exact", linestyle=:dash, lw=2)
    ylabel!("relative error")
    xlabel!("Number of colors")

    plot!(size=(400, 300), ylim=(0.95, 2.5), xlim=(0.0, 150), thickness_scaling=1.25,
        legend=:topright, left_margin=-3Plots.mm, bottom_margin=-3Plots.mm)
    savefig("colors-lp.pdf")
end

function make_flow_plot()
    plot()
    results = joinpath(pwd(), "Lifted Inference", "results", "flow.tsv")
    df = CSV.read(results, DataFrame; delim='\t')

    for x in DATASETS_FLOW
        df_x = subset(df, :name => ByRow(==(x)))
        baseline = subset(df_x, :colors => ByRow(==(0)))
        baseline = baseline[1, :]
        df_x = subset(df_x, :colors => ByRow(!=(0)))

        df_x.time = df_x.time / baseline.time * 100
        df_x.objective .= (df_x.objective ./ baseline.objective)
        df_x.objective .= max.(df_x.objective, 1 ./ df_x.objective)
        plot!(df_x.time, df_x.objective, label=splitext(x)[1], lw=2, markershape=:auto)
    end

    hline!([1.0], label="exact", linestyle=:dash, lw=2)
    ylabel!("relative error")
    xlabel!("Runtime (% of baseline)")

    plot!(size=(400, 300), ylim=(0.95, 2.0), xlim=(0.0, 1.5), thickness_scaling=1.25,
        legend=:topright, left_margin=-3Plots.mm, bottom_margin=-3Plots.mm)
    savefig("flow.pdf")
end

function make_flow_color_plot()
    plot()
    results = joinpath(pwd(), "Lifted Inference", "results", "flow.tsv")
    df = CSV.read(results, DataFrame; delim='\t')

    for x in DATASETS_FLOW
        df_x = subset(df, :name => ByRow(==(x)))
        baseline = subset(df_x, :colors => ByRow(==(0)))
        baseline = baseline[1, :]
        df_x = subset(df_x, :colors => ByRow(!=(0)))

        df_x.time = df_x.time / baseline.time * 100
        df_x.objective .= (df_x.objective ./ baseline.objective)
        df_x.objective .= max.(df_x.objective, 1 ./ df_x.objective)
        plot!(df_x.colors, df_x.objective, label=splitext(x)[1], lw=2, markershape=:auto)
    end

    hline!([1.0], label="exact", linestyle=:dash, lw=2)
    ylabel!("relative error")
    xlabel!("Number of colors")

    plot!(size=(400, 300), ylim=(0.95, 2.0), xlim=(0.0, 35), thickness_scaling=1.25,
        legend=:outerright, left_margin=-3Plots.mm, bottom_margin=-3Plots.mm)
    savefig("colors-flow.pdf")
end


function make_centrality_plot()
    plot()
    results = joinpath(pwd(), "Lifted Inference", "results", "centrality.tsv")
    df = CSV.read(results, DataFrame; delim='\t')

    for x in unique(df.name)
        df_x = subset(df, :name => ByRow(==(x)))

        plot!(df_x.time, df_x.corr, label=splitext(x)[1], lw=2, markershape=:auto)
    end

    hline!([1.0], label="exact", linestyle=:dash, lw=2)
    ylabel!("Spearman rank correlation")
    xlabel!("Runtime (% of baseline)")

    plot!(size=(400, 300), ylim=(0.80, 1.01), xlim=(0, 3), thickness_scaling=1.25,
        legend=:bottomright, left_margin=-3Plots.mm, bottom_margin=-3Plots.mm)
    savefig("centrality.pdf")
end

function make_centrality_color_plot()
    plot()
    results = joinpath(pwd(), "Lifted Inference", "results", "centrality.tsv")
    df = CSV.read(results, DataFrame; delim='\t')

    for x in unique(df.name)
        df_x = subset(df, :name => ByRow(==(x)))

        plot!(df_x.colors, df_x.corr, label=splitext(x)[1], lw=2, markershape=:auto)
    end

    hline!([1.0], label="exact", linestyle=:dash, lw=2)
    ylabel!("Spearman rank correlation")
    xlabel!("Number of colors")

    plot!(size=(400, 300), ylim=(0.80, 1.01), xlim=(0, 150), thickness_scaling=1.25,
        legend=:bottomright, left_margin=-3Plots.mm, bottom_margin=-3Plots.mm)
    savefig("colors-centrality.pdf")
end

make_lp_plot()
make_lp_color_plot()
make_flow_plot()
make_flow_color_plot()
make_centrality_plot()
make_centrality_color_plot()
