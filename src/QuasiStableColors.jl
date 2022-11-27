"""Compute quasi-stable colorings for graphs. Includes approximating max-flow, linear
programs and betweenness centrality."""
module QuasiStableColors

include("centrality.jl")
include("flow.jl")
include("optimize.jl")

export refine_stable, q_color, refine_bipartite, Color
export approx_betweenness_centrality
export lifted_maxflow
export lifted_maximize, lifted_minimize, maximize, minimize

include("misc.jl")

using Graphs
using SparseArrays
using Statistics: mean, median

Color = UInt

"""
    refine_stable(G::AbstractGraph{T})

Compute the stable coloring for the undirected graph `G`. Provided for comparasion.
"""
refine_stable(G::AbstractGraph{T}; args...) where {T} =
    q_color(G; q=0.0, args...)


"""Book-keeping datastructure for the q-coloring algorithm."""
struct ColorStats
    v::Int
    n::Int
    neighbor_base::Matrix{Float64}
    upper_base::Matrix{Float64}
    lower_base::Matrix{Float64}
    counts_base::Matrix{UInt}
    errors_base::Matrix{Float64}
end

function ColorStats(v::Int, n::Int)
    return ColorStats(
        v,
        n,
        zeros(v, n),
        zeros(n, n),
        fill(+Inf, n, n),
        zeros(n, n),
        zeros(n, n),
    )
end


"""Resize `old` to have size `n`."""
function ColorStats(old::ColorStats,v::Int, n::Int)
    result = ColorStats(
        v,
        n,
        zeros(v, n),
        zeros(n, n),
        fill(+Inf, n, n),
        zeros(n, n),
        zeros(n, n),
    )

    m = old.n
    result.neighbor_base[:, 1:m] .= old.neighbor_base
    result.upper_base[1:m, 1:m] .= old.upper_base
    result.lower_base[1:m, 1:m] .= old.lower_base
    result.counts_base[1:m, 1:m] .= old.counts_base
    result.errors_base[1:m, 1:m] .= old.errors_base

    return result
end


function pick_witness(P, P_sparse::SparseMatrixCSC{Float64,Int},
    weights::SparseMatrixCSC{Float64,Int}, stats::ColorStats)
    _, m = size(P_sparse)

    neighbor::SparseMatrixCSC{Float64,Int} = weights * P_sparse

    upper_deg = @view stats.upper_base[1:m, 1:m]
    lower_deg = @view stats.lower_base[1:m, 1:m]
    errors = @view stats.errors_base[1:m, 1:m]

    # group the rows by partition
    for i in eachindex(P)
        X::Vector{Int64} = P[i]
        upper_deg[i, :] .= transpose(maximum(neighbor[X, :], dims=1))
        lower_deg[i, :] .= transpose(minimum(neighbor[X, :], dims=1))
    end

    #errors .= (upper_deg - lower_deg) .* transpose([(length(P_i)) for P_i in P]) .* [length(P_i) for P_i in P]
    errors .= upper_deg - lower_deg

    # check for NaN or Inf
    @assert all(isfinite, errors)

    _, witness = findmax(errors)
    q_error = maximum(upper_deg - lower_deg)
    error_val = mean(errors)
    witness_i, witness_j = witness[1], witness[2]

    split_deg = mean(neighbor[P[witness_i], witness_j])

    if length(P) % 10 == 0
        @debug begin
            m = median(length(p) for p in P)
            error_sum = round(log(10, error_val), digits=2)
            "$(length(P)) colors, q-error $q_error, error sum 10^$error_sum, median size $m"
        end
    end

    return witness_i, witness_j, split_deg, error_val, q_error
end

"""
    q_color(
        G::AbstractGraph{T},
        q = 0.0,
        n_colors = Inf,
        weights::SparseMatrixCSC{<:Number,Int} = nothing,
        special =  Set{T}(),
        warm_start = Vector{Vector{T}}(),
    )


Compute a quasi-stable coloring for the undirected graph `G`. Typically, you should 
set one of:
- **`q`**: maximum q-error allowed
- **`n_colors`**: number of colors to use

Optional parameters:
- **`warm_start`**: coloring to refine. If not provided, start using trivial (single color) partitioning assumed.
- **`weights`**: edge weights to use
"""
function q_color(G::AbstractGraph{T};
    weights::Union{SparseMatrixCSC{<:Number,Int},Nothing}=nothing,
    special::Set{T}=Set{T}(),
    warm_start::Vector{Vector{T}}=Vector{Vector{T}}(),
    n_colors=Inf,
    q::Float64=0.0) where {T}

    V = Set(vertices(G))
    local P::Vector{Vector{T}}
    if length(warm_start) == 0
        P = Vector(vcat([collect(setdiff(Set(V), special))],
            [Vector([x]) for x in special]))
    else
        P = deepcopy(warm_start)
    end

    if weights === nothing
        weights::SparseMatrixCSC{Float64,Int} = copy(adjacency_matrix(G, Float64; dir=:both))
        weights.nzval .= 1.0
    end

    color_stats = ColorStats(nv(G), floor(Int, min(n_colors, 128)))

    while length(P) < n_colors
        if length(P) == color_stats.n
            color_stats = ColorStats(color_stats, nv(G), color_stats.n * 2)
        end

        P_sparse = partition_matrix(P)
        witness_i, witness_j, split_deg, _, q_error = pick_witness(P, P_sparse, weights,
            color_stats)

        if q_error <= q
            break
        end

        # split the witness_i-th color 
        retained::Vector{T} = []
        ejected::Vector{T} = []
        neighbor::SparseMatrixCSC{Float64,Int} = weights * P_sparse
        for v in P[witness_i]
            if neighbor[v, witness_j] > split_deg
                push!(ejected, v)
            else
                push!(retained, v)
            end
        end

        @assert length(retained) != 0
        @assert length(ejected) != 0

        # update the partitions
        P[witness_i] = retained
        push!(P, ejected)
    end

    P_sparse = partition_matrix(P)
    _, _, _, error, q_error = pick_witness(P, P_sparse, weights, color_stats)
    @debug "refined and got $(length(P)) colors with $q_error q-error, $error sum error"
    return P
end

"""
Equivalent to `refine_fixpoint` but optimized for bipartite graphs. Faster but less
general.
"""
function refine_bipartite(M; n_colors=Inf,
    q::Float64=0.0)
    M_t = transpose(M)
    m, n = size(M)
    P_row::Vector{Vector{Int}} = [collect(1:m-1), [m]]
    P_col::Vector{Vector{Int}} = [collect(1:n-1), [n]]

    i = 1
    err = Inf
    while i < n_colors && err > q
        if i % 10 == 0
            @debug "error, colors:" err i
        end
        row_err, row_i, col_j = _refine_bipartite(M, P_row, P_col)
        col_err, col_i, row_j = _refine_bipartite(M_t, P_col, P_row)
        old_color::Vector{Int} = []
        new_color::Vector{Int} = []
        err = max(col_err, row_err)

        if err <= q
            break
        end

        if col_err > row_err || length(P_row) == m
            # split col
            to_split = P_col[col_i]
            split_by = P_row[row_j]
            vals = sum(M_t[to_split, split_by], dims=2)
            threshold = mean(vals)

            for i in eachindex(to_split)
                if vals[i] > threshold
                    push!(new_color, to_split[i])
                else
                    push!(old_color, to_split[i])
                end
            end

            @assert length(new_color) != 0
            @assert length(old_color) != 0

            # update color
            P_col[col_i] = old_color
            push!(P_col, new_color)
        else
            # split row
            to_split = P_row[row_i]
            split_by = P_col[col_j]
            vals = sum(M[to_split, split_by], dims=2)
            threshold = mean(vals)

            for i in eachindex(to_split)
                if vals[i] > threshold
                    push!(new_color, to_split[i])
                else
                    push!(old_color, to_split[i])
                end
            end

            @assert length(new_color) != 0
            @assert length(old_color) != 0

            # update color
            P_row[row_i] = old_color
            push!(P_row, new_color)
        end
        i += 1
    end

    return P_row, P_col
end

function _refine_bipartite(M, P::Vector{Vector{Int}}, G::Vector{Vector{Int}})
    W = partition_matrix(G)
    sums = M * W
    errors = zeros(Float64, length(P), length(G))
    for i in eachindex(P)
        errors[i, :] = (maximum(sums[P[i], :], dims=1) - minimum(sums[P[i], :], dims=1)) * length(P[i])
    end
    val, witness = findmax(errors)
    witness_i, witness_j = witness[1], witness[2]
    return val, witness_i, witness_j
end

function geomean(a)
    s = 0.0
    n = length(a)
    for i = 1:n
        @inbounds s += log(a[i])
    end
    return exp(s / n)
end
end
