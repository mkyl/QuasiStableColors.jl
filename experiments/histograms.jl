using Plots

include("../main.jl")
include("../optimize.jl")
include("../datasets.jl")
include("measure_lp.jl")

function dblp()
    G = datasets.dblp()

    P = refine_fixpoint(G, early_stop=200)
    L = [length(x) for x in P]
    L_log = log10.(L)

    histogram(L_log, bins=5, xlabel="Partition size 10ˣ", ylabel="Count", legend=false,
        size=(400, 125), bottom_margin=3Plots.mm)
    savefig("histogram-dblp.pdf")
end

function epinions()
    G = datasets.epinions()
    P = refine_fixpoint(G, early_stop=100)
    L = [length(x) for x in P]
    L_log = log10.(L)

    histogram(L_log, bins=5, xlabel="Partition size 10ˣ", ylabel="Count", legend=false,
        size=(400, 125), bottom_margin=3Plots.mm)
    savefig("histogram-epinions.pdf")
end

function cells()
    G, C, _, _ = datasets.cells()
    P = refine_fixpoint(G, weights=C, early_stop=100)
    L = [length(x) for x in P]
    L_log = log10.(L)

    histogram(L_log, bins=5, xlabel="Partition size 10ˣ", ylabel="Count", legend=false,
        size=(400, 125), bottom_margin=3Plots.mm)
    savefig("histogram-cells.pdf")
end

function nug()
    A, b, c = mps_to_matrices("Datasets/optimization/nug08-3rd.mps")
    Ã::SparseMatrixCSC{Float64,Int} = [A b; transpose(c) 0]
    G, W, _, _ = matrix_to_graph(Ã)
    P = refine_fixpoint(G, weights=W, early_stop=100)
    L = [length(x) for x in P]
    L_log = log10.(L)

    histogram(L_log, bins=5, xlabel="Partition size 10ˣ", ylabel="Count", legend=false,
        size=(400, 125), bottom_margin=3Plots.mm)
    savefig("histogram-nug08-3rd.pdf")
end

nug()
