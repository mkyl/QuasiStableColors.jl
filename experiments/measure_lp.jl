using QPSReader
using SparseArrays
using CSV
using DataFrames

DATASETS = [
    "supportcase10.mps",
    "nug08-3rd.mps",
    "qap15.mps",
    "ex10.mps",
]

"""Read MPS file and return A, b, c from the system `min c^T x` where `A x ≥ b`."""
function mps_to_matrices(filename)
    mps = readqps(filename)
    I, J, V = mps.arows, mps.acols, mps.avals
    c = mps.c
    b_lo = mps.lcon
    A = sparse(I, J, V, length(b_lo), length(c))

    # drop constraint rows of the format ≥-∞
    no_op = findall(!=(-Inf), b_lo)
    A_lo = A[no_op, :]
    b_lo = b_lo[no_op]

    # upper bounds
    b_up = mps.ucon
    no_op = findall(!=(Inf), b_up)
    A_up = A[no_op, :]
    b_up = b_up[no_op]

    # convert to lower bounds
    A_up = -1 * A_up
    b_up = -1 * b_up

    # combine constraints
    A = vcat(A_lo, A_up)
    b = vcat(b_lo, b_up)

    return A, b, c
end

include("../optimize.jl")
ENV["JULIA_DEBUG"] = Main

function measure()
    parent = pwd()
    mkpath(joinpath(parent, "results"))
    results = joinpath(parent, "Lifted Inference", "results", "lp.tsv")
    if isfile(results)
        df = CSV.read(results, DataFrame; delim='\t')
    else
        df = DataFrame(name=String[], colors=Int[], time=Float64[], size=Int[],
            objective=Float64[])
    end

    for x in DATASETS
        A, b, c = mps_to_matrices("Datasets/optimization/" * x)

        # for the JIT
        try
            lifted_minimize(A, b, c; early_stop=10)
        catch
        end

        for i in [4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200]
            #for i in 70:120
            if !any((df.name .== x) .& (df.colors .== i))
                try
                    stats = @timed lifted_minimize(A, b, c; early_stop=i)
                    z = stats.value
                    if z !== nothing
                        time = stats.time
                        size = 1
                        obj = z
                        push!(df, [x i time size obj])
                        CSV.write(results, df; delim='\t')
                    end
                catch LoadError
                end
            end
        end

        if !any((df.name .== x) .& (df.colors .== 0))
            stats = @timed minimize(A, b, c)
            z_star = stats.value
            size = 1
            time = stats.time
            obj = transpose(c) * z_star
            push!(df, [x 0 time size obj])
            CSV.write(results, df; delim='\t')
        end
    end
end

measure()
