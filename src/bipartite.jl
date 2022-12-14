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
