struct KuramotoSys{T1<:AbstractSparseMatrix{<:Real},
                   T2<:AbstractVector{<:Real}} <: ODESys
    A::T1
    ω::T2

    function KuramotoSys(A::T1, ω::T2) where {T1 <: AbstractSparseMatrix{<:Real},
                                              T2 <: AbstractVector{<:Real}}
        nr, nc = size(A)
        nr == nc || error("Interaction matrix a must be square")

        length(ω) == nr || error("Dimensions of ω must be compatible with A")
        new{T1, T2}(A, ω)
    end

end

# arbitrary (not-necessarily sparse) input for A
KuramotoSys(A::AbstractMatrix, ω) = KuramotoSys(sparse(A), ω)

# zero-frequency case
KuramotoSys(A::AbstractMatrix{T}) where {T} = 
    KuramotoSys(A, zeros(T, size(A, 2)))

# from a graph
KuramotoSys(G::AbstractGraph, ω) = KuramotoSys(adjacency_matrix(G), ω)

# zero-frequency case
KuramotoSys(G::AbstractGraph) = KuramotoSys(adjacency_matrix(G))

state_size(sys::KuramotoSys) = size(sys.ω)

function rhs(dθdt, θ, sys::KuramotoSys, t=0.0)
    A = sys.A
    rows = rowvals(A)
    nzv = nonzeros(A)

    copy!(dθdt, sys.ω)
    @inbounds for j in 1:size(A, 2), r in nzrange(A, j)
        i = rows[r]
        dθdt[j] += nzv[r]*sin(θ[i] - θ[j])
    end
end

function jac(J::AbstractMatrix{T}, θ, sys::KuramotoSys, t=0.0) where {T <: Real}
    A = sys.A
    rows = rowvals(A)
    nzv = nonzeros(A)

    fill!(J, zero(T))
    @inbounds for j in 1:size(A, 2), r in nzrange(A, j)
        i = rows[r]
        i != j || continue
        v = nzv[r]*cos(θ[i] - θ[j])
        J[j, i] = v
        J[j, j] -= v
    end
end
