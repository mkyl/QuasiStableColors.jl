include("main.jl")

function main()
    V, E = dblp()
    # println(length(refine_stable(V, E)))
    println(length(refine_abs_approx(V, E, 2.0)))
    # println(length(refine_rel_approx(V, E, 2.0)))
end

function plot_robustness!()
    interval = 1:5
    M = zeros(Int, (length(interval), length(Datasets) * 3))

    for (i, D) in enumerate(Datasets)
        V, E = D()
        col_abs = map(x -> length(refine_abs_approx(V, E, 2^x)), interval)
        col_rel = map(x -> length(refine_rel_approx(V, E, 2^x)), interval)
        col_stab = fill(length(refine_stable(V, E)), length(interval))
        M[:, i] = col_abs
        M[:, length(Datasets)+i] = col_rel
        M[:, 2*length(Datasets)+i] = col_stab
    end
    plot(map(x -> 2^x, interval), M, layout = (length(Datasets), 1),
        legend = [true false],
        seriestype = :path,
        yscale = :log10, marker = true, title = ["OpenFlights" "DBLP"],
        titlelocation = :left, xlabel = ["" "Threshold (ε or c)"],
        ylabel = ["" ".                              Number of colors"],
        labels = ["±ε-stable" "stable" "÷c-stable"], width = 5, size = (400, 400),
        markersize = 6)
    savefig("Lifted Inference/paper/graphics/robustness.pdf")
    return M
end
