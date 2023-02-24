include("misc.jl")

using Graphs
using SparseArrays
using Statistics: mean, median

Color = UInt

const BASE_MATRIX_SIZE = 128

"""
    refine_stable(G::AbstractGraph{T})

Compute the stable coloring for the undirected graph `G`. Provided for comparasion.
"""
refine_stable(G::AbstractGraph{T}; args...) where {T} =
    q_color(G; q=0.0, args...)


"""Book-keeping datastructure for the q-coloring algorithm."""
mutable struct ColorStats
    v::Int
    n::Int
    neighbor::SparseMatrixCSC{Float64,Int}
    upper_base::Matrix{Float64}
    lower_base::Matrix{Float64}
    counts_base::Matrix{UInt}
    errors_base::Matrix{Float64}
end

"""Output of a quasi-stable coloring algorithm."""
struct QuasiStableColoring{T<:Int}
    partitions::Vector{Vector{T}}
    mapping::Dict{T,Color}
    max_q::Number
    avg_q::Number
end

function ColorStats(v::Int, n::Int)
    return ColorStats(
        v,
        n,
        spzeros(v, n),
        zeros(n, n),
        fill(+Inf, n, n),
        zeros(n, n),
        zeros(n, n),
    )
end


"""Resize `old` to have size `n`."""
function ColorStats(old::ColorStats, v::Int, n::Int)
    result = ColorStats(
        v,
        n,
        old.neighbor,
        zeros(n, n),
        fill(+Inf, n, n),
        zeros(n, n),
        zeros(n, n),
    )

    m = old.n
    result.upper_base[1:m, 1:m] .= old.upper_base
    result.lower_base[1:m, 1:m] .= old.lower_base
    result.counts_base[1:m, 1:m] .= old.counts_base
    result.errors_base[1:m, 1:m] .= old.errors_base

    return result
end


"""Update stats by recomputing everything."""
function update_stats!(stats::ColorStats, weights::SparseMatrixCSC{<:Number,Int},
    P::Vector{Vector{T}}; weighting=false) where {T}
    P_sparse = partition_matrix(P)
    stats.neighbor = weights * P_sparse

    m = length(P)
    upper_deg = @view stats.upper_base[1:m, 1:m]
    lower_deg = @view stats.lower_base[1:m, 1:m]
    errors = @view stats.errors_base[1:m, 1:m]

    # group the rows by partition
    for i in eachindex(P)
        X::Vector{Int64} = P[i]
        upper_deg[i, :] .= transpose(maximum(stats.neighbor[X, :], dims=1))
        lower_deg[i, :] .= transpose(minimum(stats.neighbor[X, :], dims=1))
    end

    if weighting
        errors .= (upper_deg - lower_deg) .* transpose([(length(P_i)) for P_i in P])
        # .* [length(P_i) for P_i in P]
    else
        errors .= upper_deg - lower_deg
    end

    # check for NaN or Inf
    @assert all(isfinite, errors)
end

"""Update stats by splitting one partition."""
function update_stats!(stats::ColorStats, weights::SparseMatrixCSC{<:Number,Int},
    P::Vector{Vector{T}}, old::Int, new::Int; weighting=weighting) where {T}
    old_nodes = P[old]
    new_nodes = P[new]

    # expand neighbor by one column
    stats.neighbor = [stats.neighbor spzeros(stats.v, 1)]

    # update columns for old and new color
    old_degs = sum(weights[:, old_nodes], dims=2)
    new_degs = sum(weights[:, new_nodes], dims=2)
    stats.neighbor[:, old] = old_degs
    stats.neighbor[:, new] = new_degs

    m = length(P)
    upper_deg = @view stats.upper_base[1:m, 1:m]
    lower_deg = @view stats.lower_base[1:m, 1:m]
    errors = @view stats.errors_base[1:m, 1:m]

    upper_deg[old, :] .= transpose(maximum(stats.neighbor[P[old], :], dims=1))
    lower_deg[old, :] .= transpose(minimum(stats.neighbor[P[old], :], dims=1))

    upper_deg[new, :] .= transpose(maximum(stats.neighbor[P[new], :], dims=1))
    lower_deg[new, :] .= transpose(minimum(stats.neighbor[P[new], :], dims=1))

    for i in eachindex(P)
        X::Vector{Int64} = P[i]
        upper_deg[i, old] = maximum(stats.neighbor[X, old])
        lower_deg[i, old] = minimum(stats.neighbor[X, old])

        upper_deg[i, new] = maximum(stats.neighbor[X, new])
        lower_deg[i, new] = minimum(stats.neighbor[X, new])
    end

    if weighting
        errors .= (upper_deg - lower_deg) .* transpose([(length(P_i)) for P_i in P])
        # .* [length(P_i) for P_i in P]
    else
        errors .= upper_deg - lower_deg
    end

    # check for NaN or Inf
    @assert all(isfinite, errors)
end

"""Split the color `U` according degree into the color `V`."""
function split_color!(P::Vector{Vector{T}}, stats::ColorStats, U::Int,
    V::Int, threshold::Float64) where {T}
    retained::Vector{T} = []
    ejected::Vector{T} = []
    for x in P[U]
        if stats.neighbor[x, V] > threshold
            push!(ejected, x)
        else
            push!(retained, x)
        end
    end

    @assert length(retained) != 0
    @assert length(ejected) != 0

    # update the partitions
    P[U] = retained
    push!(P, ejected)
end

function pick_witness(P, stats::ColorStats)
    m = length(P)
    upper_deg = @view stats.upper_base[1:m, 1:m]
    lower_deg = @view stats.lower_base[1:m, 1:m]
    errors = @view stats.errors_base[1:m, 1:m]

    _, witness = findmax(errors)
    q_error = maximum(upper_deg - lower_deg)
    witness_i, witness_j = witness[1], witness[2]

    split_deg = mean(stats.neighbor[P[witness_i], witness_j])

    if length(P) % 10 == 0
        @debug begin
            m = median(length(p) for p in P)
            error_sum = round(log(10, error_val), digits=2)
            "$(length(P)) colors: max $q_error err, err sum 10^$error_sum, median $m"
        end
    end

    return witness_i, witness_j, split_deg, q_error
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


Compute a quasi-stable coloring for the graph `G`. Typically, you should 
set one of:
- **`q`**: maximum q-error allowed
- **`n_colors`**: number of colors to use

Advanced, optional parameters:
- **`warm_start`**: coloring to refine. If not provided, start using trivial 
(single color) partitioning assumed.
- **`weights`**: edge weights to use
"""
function q_color(G::AbstractGraph{T};
    weights::Union{SparseMatrixCSC{<:Number,Int},Nothing}=nothing,
    special::Set{T}=Set{T}(),
    warm_start::Vector{Vector{T}}=Vector{Vector{T}}(),
    weighting=false,
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
        weights::SparseMatrixCSC{Float64,Int} = copy(adjacency_matrix(G, Float64))
        weights.nzval .= 1.0
    end
    weightsᵀ = SparseMatrixCSC(transpose(weights))

    # + symbol represents out-degree, - symbol in-degree 
    color_stats⁺ = ColorStats(nv(G), floor(Int, min(n_colors, BASE_MATRIX_SIZE)))
    color_stats⁻ = ColorStats(nv(G), floor(Int, min(n_colors, BASE_MATRIX_SIZE)))
    update_stats!(color_stats⁺, weights, P, weighting=weighting)
    update_stats!(color_stats⁻, weightsᵀ, P, weighting=weighting)

    while length(P) < n_colors
        # check if we need to grow data structures
        if length(P) == color_stats⁺.n
            color_stats⁺ = ColorStats(color_stats⁺, nv(G), color_stats⁺.n * 2)
            color_stats⁻ = ColorStats(color_stats⁻, nv(G), color_stats⁻.n * 2)
        end

        witness⁺ᵢ, witness⁺ⱼ, split_deg⁺, q_error⁺ = pick_witness(P, color_stats⁺)
        witness⁻ᵢ, witness⁻ⱼ, split_deg⁻, q_error⁻ = pick_witness(P, color_stats⁻)

        if q_error⁺ ≤ q && q_error⁻ ≤ q
            break
        end

        if q_error⁺ ≥ q_error⁻
            split_color!(P, color_stats⁺, witness⁺ᵢ, witness⁺ⱼ, split_deg⁺)
            update_stats!(color_stats⁺, weights, P, witness⁺ᵢ, length(P), weighting=weighting)
            update_stats!(color_stats⁻, weightsᵀ, P, witness⁺ᵢ, length(P), weighting=weighting)
        else
            split_color!(P, color_stats⁻, witness⁻ᵢ, witness⁻ⱼ, split_deg⁻)
            update_stats!(color_stats⁺, weights, P, witness⁻ᵢ, length(P), weighting=weighting)
            update_stats!(color_stats⁻, weightsᵀ, P, witness⁻ᵢ, length(P), weighting=weighting)
        end
    end

    _, _, _, q_error⁺ = pick_witness(P, color_stats⁺)
    _, _, _, q_error⁻ = pick_witness(P, color_stats⁻)
    q_error = max(q_error⁺, q_error⁻)
    @debug "refined and got $(length(P)) colors with $q_error q-error"
    return P
end
