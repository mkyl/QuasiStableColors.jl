"""Functions related to linear optimization."""
module Optimize

export lifted_maximize, lifted_minimize, maximize, minimize

using Graphs
using JuMP, Tulip
using MathOptInterface
using SparseArrays
using LinearAlgebra

using QuasiStableColors

"""Approximate linear program `min c^T x where A x >= b, x >=0`."""
lifted_minimize(A, b, c; args...) = _lifted_opt(A, b, c; obj=MathOptInterface.MIN_SENSE,
    args...)
"""Approximate linear program `max c^T x where A x <= b, x>=0`."""
lifted_maximize(A, b, c; args...) = _lifted_opt(A, b, c; obj=MathOptInterface.MAX_SENSE,
    args...)

function _lifted_opt(A, b::Vector, c::Vector; obj=MathOptInterface.MIN_SENSE, args...)
    Ã::SparseMatrixCSC{Float64,Int} = [A b; transpose(c) 0]
    m, n = size(Ã)

    P_row, P_col = refine_bipartite(Ã; args...)
    P_row = filter!(x -> x ≠ [m], P_row)
    P_col = filter!(x -> x ≠ [n], P_col)

    U, V = P_to_UV(P_row), P_to_UV(P_col)

    @assert V * transpose(V) ≈ I
    A₂ = U * A * transpose(V)
    b₂ = U * b
    c₂ = V * c

    @info "lifted linear program size: " size(A₂)
    # solve the lifted problem
    z₂ = _optimize(A₂, b₂, c₂; obj=obj)
    if (z₂ !== nothing)
        c₂' * z₂
    else
        nothing
    end
end

"Minimize c^T x where A x >= b, x>=0."
minimize(A, b, c) = _optimize(A, b, c; obj=MathOptInterface.MIN_SENSE)

"Maximize c^T x where A x <= b, x>=0."
maximize(A, b, c) = _optimize(A, b, c; obj=MathOptInterface.MAX_SENSE)

function _optimize(A, b, c; obj=MathOptInterface.MAX_SENSE)
    @assert typeof(b) <: Vector
    @assert typeof(c) <: Vector
    @assert size(A) == (length(b), length(c))

    model = Model(Tulip.Optimizer)
    @variable(model, x[1:length(c)] >= 0)
    @objective(model, obj, c' * x)
    if obj == MathOptInterface.MAX_SENSE
        @constraint(model, A * x .<= b)
    else
        @constraint(model, A * x .>= b)
    end
    optimize!(model)
    if termination_status(model) == MOI.OPTIMAL
        return value.(x)
    else
        @warn "Solver did not find global minimum" termination_status(model)
        return nothing
    end
end

"""
Convert matrix M to a bipartite graph G=(U, V, E), where U, V are rows, columns
respectively.
"""
function matrix_to_graph(M)
    n, m = size(M)
    G = SimpleDiGraph(n + m)
    I, J, V = [], [], []

    I_m, J_m, V_m = findnz(M)
    for x in 1:length(V_m)
        i = I_m[x]
        j = J_m[x]
        v = V_m[x]
        add_edge!(G, i, n + j)

        # same as W[i, n+j] = M[i, j]
        push!(I, i)
        push!(J, n + j)
        push!(V, v)

        # same as W[n+j, i] = M[i, j]
        push!(I, n + j)
        push!(J, i)
        push!(V, v)
    end
    W::SparseMatrixCSC{Float64,Int64} = sparse(I, J, V, n + m, n + m)
    U = Set(range(1, n))
    V = Set(range(n + 1, n + m))
    @assert is_bipartite(G)
    return G, W, U, V
end

"""Transform partition P into the conversion matrix U or V"""
function P_to_UV(P)
    n = length(P)
    m = sum(length(x) for x in P)
    U::Matrix = zeros(n, m)
    for (i, X) in enumerate(P)
        for x in X
            U[i, x] = 1 / sqrt(length(X))
        end
    end
    U
end

end