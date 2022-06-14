# Utility Functions
function partition_matrix(P::Vector{Vector{T}}) where {T}
    n = sum(length(X) for X in P)
    i = 1
    I::Vector{Int}, J::Vector{Int}, V::Vector{Float64} = zeros(n), zeros(n), zeros(n)
    for (c, X) in enumerate(P)
        for x in X
            I[i] = x
            J[i] = c
            V[i] = 1.0
            i += 1
        end
    end

    P_sparse::SparseMatrixCSC{Float64,Int} = sparse(I, J, V, n, length(P))
    return P_sparse
end