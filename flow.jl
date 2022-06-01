include("main.jl")

using Graphs, GraphsFlows, StatsBase

function uplo(x::AbstractVector; prop::Real=0.0, count::Integer=0)
    n = length(x)
    n > 0 || throw(ArgumentError("x can not be empty."))

    if count == 0
        0 <= prop < 0.5 || throw(ArgumentError("prop must satisfy 0 ≤ prop < 0.5."))
        count = floor(Int, n * prop)
    else
        prop == 0 || throw(ArgumentError("prop and count can not both be > 0."))
        0 <= count < n / 2 || throw(ArgumentError("count must satisfy 0 ≤ count < length(x)/2."))
    end

    # indices for lowest count values
    x2 = copy(x)
    #lo = partialsort!(x2, 1:count+1)[end]
    lo = 0
    # indices for largest count values
    up = partialsort!(x2, n-count:n)[1]

    up, lo
end

# from https://github.com/JuliaStats/StatsBase.jl/
# blob/08d4b77a4b42ef8cadb67b98da65e7fbd9959e0b/src/robust.jl#L32-L52
function trim(x::AbstractVector; prop::Real=0.0, count::Integer=0)
    up, lo = uplo(x; prop=prop, count=count)
    (xi for xi in x if lo <= xi <= up)
end

function quotient_graph(G::AbstractGraph{T};
    weights::Union{SparseMatrixCSC{<:Number,Int},Nothing}=nothing, args...) where {T}
    P = refine_fixpoint(G; weights=weights, args...)

    k = length(P)
    G₀ = SimpleDiGraph(k)

    𝜑 = Dict{T,Color}()
    sizehint!(𝜑, nv(G))
    for (color, nodes) in enumerate(P)
        for x in nodes
            𝜑[x] = color
        end
    end

    P_sparse = partition_matrix(P)
    neighbor = weights * P_sparse

    @assert sum(neighbor, dims=2) == sum(weights, dims=2)

    C₀::Matrix{Float64} = zeros(k, k)
    for i in eachindex(P)
        X::Vector{Int64} = P[i]
        C₀[i, :] = sum(neighbor[X, :], dims=1)
        #C₀[i, :] = [mean(trim(neighbor[X, j], prop=0.3)) for j in eachindex(P)] * length(X)
    end

    # @assert sum(C₀) <= sum(weights)

    for i in findall(w -> w != 0.0, C₀)
        u, v = i[1], i[2]
        add_edge!(G₀, u, v)
    end

    return G₀, C₀, 𝜑, P
end

function lifted_maxflow(G, s::Int, t::Int; args...)::Number
    G₀, C₀, 𝜑, _ = quotient_graph(G; special=Set([s, t]), args...)
    s₀, t₀ = 𝜑[s], 𝜑[t]
    f₀, _ = maximum_flow(G₀, s₀, t₀, C₀)
    return f₀
end
