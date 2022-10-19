"""Compute quasi-stable colorings for graphs. Includes approximating max-flow, linear
programs and betweenness centrality."""
module QuasiStableColors

include("centrality.jl")
include("flow.jl")
include("optimize.jl")

export refine_stable, refine_fixpoint, refine_bipartite, Color
export approx_betweenness_centrality
export lifted_maxflow
export lifted_maximize, lifted_minimize, maximize, minimize

include("misc.jl")

using Graphs
using SparseArrays
using DataStructures: counter
using Statistics: mean, median

Color = UInt

"""
Compute the stable coloring for the graph with vertices `V` and edges `E`.
"""
refine_stable(G::AbstractGraph{T}; args...) where {T} =
    refine_fixpoint(G; eps=0.0, args...)


function pick_witness(P, P_sparse, weights, neighbor_base, upper_base, lower_base,
    counts_base, errors_base)
    n, m = size(P_sparse)

    neighbor::SparseMatrixCSC{Float64,Int} = weights * P_sparse

    upper_deg = @view upper_base[1:m, 1:m]
    lower_deg = @view lower_base[1:m, 1:m]
    errors = @view errors_base[1:m, 1:m]
    counts = @view counts_base[1:m, 1:m]

    # group the rows by partition
    for i in eachindex(P)
        X::Vector{Int64} = P[i]
        upper_deg[i, :] .= transpose(maximum(neighbor[X, :], dims=1))
        lower_deg[i, :] .= transpose(minimum(neighbor[X, :], dims=1))
    end

    errors .= (upper_deg - lower_deg) .* transpose([(length(P_i)) for P_i in P]) .* [length(P_i) for P_i in P]
    #errors .= upper_deg - lower_deg

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
Compute a quasi-stable coloring for the graph with vertices `V` and edges `E`. Provide
`eps` for a coloring with q-error eps or `early_stop` for a coloring with a specific size.
"""
function refine_fixpoint(G::AbstractGraph{T};
    weights::Union{SparseMatrixCSC{<:Number,Int},Nothing}=nothing,
    special::Set{T}=Set{T}(),
    warm_start::Vector{Vector{T}}=Vector{Vector{T}}(),
    early_stop=Inf,
    eps::Float64=0.0) where {T}

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
        weights.nzval .= 1
    end

    if early_stop == Inf
        n = 250
    else
        n = early_stop
    end
    neighbor_base::Matrix{Float64} = zeros(nv(G), n)
    upper_base::Matrix{Float64} = zeros(n, n)
    lower_base::Matrix{Float64} = fill(+Inf, n, n)
    counts_base::Matrix{UInt} = zeros(n, n)
    errors_base::Matrix{Float64} = zeros(n, n)

    while length(P) < early_stop
        P_sparse = partition_matrix(P)
        witness_i, witness_j, split_deg, error, q_error = pick_witness(P, P_sparse, weights,
            neighbor_base, upper_base, lower_base, counts_base, errors_base)

        if q_error <= eps
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
    _, _, _, error, q_error = pick_witness(P, P_sparse, weights,
        neighbor_base, upper_base, lower_base, counts_base, errors_base)
    @debug "refined and got $(length(P)) colors with $q_error q-error, $error sum error"
    return P
end

"""
Equivalent to `refine_fixpoint()` but optimized for bipartite graphs.
"""
function refine_bipartite(M; early_stop=Inf,
    eps::Float64=0.0)
    M_t = transpose(M)
    m, n = size(M)
    P_row::Vector{Vector{Int}} = [collect(1:m-1), [m]]
    P_col::Vector{Vector{Int}} = [collect(1:n-1), [n]]

    i = 1
    err = Inf
    while i < early_stop && err > eps
        if i % 10 == 0
            @debug "error, colors:" err i
        end
        row_err, row_i, col_j = _refine_bipartite(M, P_row, P_col)
        col_err, col_i, row_j = _refine_bipartite(M_t, P_col, P_row)
        old_color::Vector{Int} = []
        new_color::Vector{Int} = []
        err = max(col_err, row_err)

        if err <= eps
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
