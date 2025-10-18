function gen_kuramoto(n, k_avg, p_neg; is_directed=false,
                      force_connected=true,
                      force_unstable=true,
                      unstable_tol=1e-3,
                      max_tries=1000)
    #p = is_directed ? k_avg/(n-1)/2 : k_avg/(n-1)
    p = k_avg / (n - 1)
    θ₀ = zeros(n)
    tries = 0

    while tries < max_tries
        tries += 1
        g = erdos_renyi(n, p; is_directed=is_directed)

        if force_connected && ~is_weakly_connected(g)
            continue
        end

        # get upper triangle of adjacency matrix (so there is only one nonzero
        # entry per edge)
        #A = convert(SparseMatrixCSC{Float64}, triu(adjacency_matrix(G)))
        A = adjacency_matrix(g, Float64)

        # change weights to uniform random numbers
        # make the weight negative with probability p_neg
        # rows = rowvals(A)
        # for j=1:n, row_idx=nzrange(A, j)
        #     i = rows[row_idx]
        #     w = rand()
        #     A[i, j] = rand() < p_neg ? -w : w
        # end
        for e in edges(g)
            i = src(e)
            j = dst(e)
            w = rand()
            w = rand() < p_neg ? -w : w
            A[i, j] = w
            if !is_directed
                A[j, i] = w
            end
        end

        # make symmetric
        #@. A = A + A'
        # A = Symmetric(A)

        sys = KuramotoSys(A)
        if force_unstable && all(λ -> real(λ) < unstable_tol, eigvals(jac(θ₀, sys)))
            continue
        end
        return sys
    end
    error("Exceeded max_tries.")
end
