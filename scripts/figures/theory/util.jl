using Distributions
using LinearAlgebra
using LoopVectorization
using Graphs
using Random
using TemporalSync
using SparseArrays

struct Directed end

struct Undirected end

function all_works(snapshots, x0)
    T = eltype(x0)
    n = length(x0)

    W = zeros(T, length(snapshots))
    
    F = zeros(T, n)
    force(KuramotoEnergy(), F, x0)

    dxdt = zeros(T, n)

    sins = zeros(T, n, n)
    for i=1:n, j=i+1:n
        v = sin(x0[i] - x0[j])
        sins[i, j] = v
        sins[j, i] = -v
    end

    for (i, s) in enumerate(snapshots)
        fill!(dxdt, zero(T))
        A = s.A
        rows = rowvals(A)
        nzv = nonzeros(A)
        for j in 1:size(A, 2), r in nzrange(A, j)
            k = rows[r]
            dxdt[j] += nzv[r]*sins[k, j]
        end

        W[i] = dot(dxdt, F)
    end
    return W
end

function weight_mean(n, k_avg, ρ)
    # mean of OUR weight distribution (incl. zero weights), which is:
    # 1) undirected links (no self-loops) created independently w/ avg. degree k_avg
    # 2) non-zero link weights assigned to be U(-1, 0) w/ prob. ρ, or
    #    from U(0, 1) w/ prob 1-ρ.
    return k_avg / (n-1)  * (1/2 - ρ)
end

function weight_var(n, k_avg, ρ)
    # variance of OUR weight distribution (incl. zero weights) as described
    # in weight_mean above.
    return k_avg / (n-1) / 3 - weight_mean(n, k_avg, ρ)^2
end

weight_std(n, k_avg, ρ) = sqrt(weight_var(n, k_avg, ρ))

function G(x, i, j)
    result = zero(eltype(x))
    n = length(x)
    for k in 1:n
        result += sin(x[k] - x[i])
    end
    result *= -2/n^2 * sin(x[j] - x[i])
    return result
end

function work_mean_and_std(mode::Directed, x0, n, k_avg, ρ)
    μ = 0.0
    σ² = 0.0
    for i=1:n, j=1:n
        i == j && continue
        tmp = G(x0, i, j)
        μ += tmp
        σ² +=  tmp^2
    end

    μ *= weight_mean(n, k_avg, ρ)
    σ² *= weight_var(n, k_avg, ρ)
    σ = sqrt(σ²)

    return μ, σ
end

function work_mean_and_std(mode::Undirected, x0, n, k_avg, ρ)
    μ = 0.0
    σ² = 0.0
    for i=1:n, j=i+1:n
        tmp = G(x0, i, j) + G(x0, j, i)
        μ += tmp
        σ² +=  tmp^2
    end

    μ *= weight_mean(n, k_avg/2, ρ)
    σ² *= weight_var(n, k_avg/2, ρ)
    σ = sqrt(σ²)

    return μ, σ
end