"""Functions related to computing flows on network graphs."""
module Flow

export lifted_maxflow

using QuasiStableColors
using Graphs, GraphsFlows, SparseArrays

include("misc.jl")

function _quotient_graph(G::AbstractGraph{T};
    weights::Union{SparseMatrixCSC{<:Number,Int},Nothing}=nothing, args...) where {T}
    P = refine_fixpoint(G; weights=weights, args...)

    k = length(P)
    Gâ = SimpleDiGraph(k)

    ğ = Dict{T,Color}()
    sizehint!(ğ, nv(G))
    for (color, nodes) in enumerate(P)
        for x in nodes
            ğ[x] = color
        end
    end

    P_sparse = partition_matrix(P)

    if weights === nothing
        @info "Assuming unit capacities"
        weights::SparseMatrixCSC{Float64,Int} = adjacency_matrix(G)
    end

    neighbor = weights * P_sparse

    @assert sum(neighbor, dims=2) == sum(weights, dims=2)

    Câ::Matrix{Float64} = zeros(k, k)
    for i in eachindex(P)
        X::Vector{Int64} = P[i]
        Câ[i, :] = sum(neighbor[X, :], dims=1)
    end

    for i in findall(w -> w != 0.0, Câ)
        u, v = i[1], i[2]
        add_edge!(Gâ, u, v)
    end

    return Gâ, Câ, ğ, P
end

"""Compute the maximum flow from `s` to `t` in flow network `G`. Capacities defined
 using `weights=`."""
function lifted_maxflow(G, s::Int, t::Int; args...)::Number
    Gâ, Câ, ğ, _ = _quotient_graph(G; special=Set([s, t]), args...)
    sâ, tâ = ğ[s], ğ[t]
    fâ, _ = maximum_flow(Gâ, sâ, tâ, Câ)
    return fâ
end
end
