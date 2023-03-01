"""Functions related to linear optimization."""
module Optimize

export lifted_maximize, lifted_minimize, maximize, minimize

using Graphs
using SparseArrays
using LinearAlgebra
using MathOptInterface
const MOI = MathOptInterface

using QuasiStableColors

"""
    lifted_minimize(
        A,
        b::Vector,
        c::Vector,
        q=0.0,
        n_colors=Inf,
    )

Approximate the linear program ``\\min c^T x \\text{ where } A x \\geq b, x \\geq 0``. 

Uses a quasi-stable coloring with maximum error `q` or `n_colors` colors, whichever
is smaller."""
lifted_minimize(optimizer, A, b, c; args...) = _lifted_opt(optimizer, A, b, c; obj=MathOptInterface.MIN_SENSE,
    args...)
"""Same as `lifted_minimize` but for the linear program
``\\max c^T x \\text{ where } A x \\leq b, x \\geq 0``."""
lifted_maximize(optimizer, A, b, c; args...) = _lifted_opt(optimizer, A, b, c; obj=MathOptInterface.MAX_SENSE,
    args...)

function _lifted_opt(optimizer::T, A, b::W, c::W; obj=MathOptInterface.MAX_SENSE, args...) where {
    T<:MOI.AbstractOptimizer,W<:Vector{Float64}}
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

    @debug "lifted linear program size: " size(A₂)
    # solve the lifted problem
    z₂ = _optimize(optimizer, A₂, b₂, c₂; obj=obj)
    if (z₂ !== nothing)
        c₂' * z₂
    else
        nothing
    end
end

"Minimize c^T x where A x >= b, x>=0."
minimize(optimizer, A, b, c) = _optimize(optimizer, A, b, c; obj=MathOptInterface.MIN_SENSE)

"Maximize c^T x where A x <= b, x>=0."
maximize(optimizer, A, b, c) = _optimize(optimizer, A, b, c; obj=MathOptInterface.MAX_SENSE)

function _optimize(optimizer::T, A::V, b::U, c::U; obj=MathOptInterface.MAX_SENSE) where {
    T<:MOI.AbstractOptimizer,U<:AbstractVector{Float64},V<:AbstractMatrix{Float64}}
    @assert size(A) == (length(b), length(c))

    x = MOI.add_variables(optimizer, length(c))
    MOI.set(
        optimizer,
        MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
        MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(c, x), 0.0),
    )

    MOI.set(optimizer, MOI.ObjectiveSense(), obj)

    # matrix of constraints
    for (i, row) in enumerate(eachrow(A))
        row_function = MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(row, x), 0.0)
        if obj == MathOptInterface.MAX_SENSE
            MOI.add_constraint(optimizer, row_function, MOI.LessThan(b[i]))
        else
            MOI.add_constraint(optimizer, row_function, MOI.GreaterThan(b[i]))
        end
    end

    # positive constraint
    for x_i in x
        MOI.add_constraint(optimizer, x_i, MOI.GreaterThan(0.0))
    end

    MOI.optimize!(optimizer)

    status = MOI.get(optimizer, MOI.TerminationStatus())
    if status == MOI.OPTIMAL
        return MOI.get(optimizer, MOI.VariablePrimal(), x)
    else
        @warn "Solver did not find global minimum" status
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