# new energy functions should implement energy and force
abstract type EnergyFunc end;
abstract type SmoothEnergyFunc <: EnergyFunc end;
abstract type NonSmoothEnergyFunc <: EnergyFunc end;

# don't broadcast over EnergyFunc
broadcastable(e::EnergyFunc) = Ref(e)

# dummy functions defining interface
function energy(::EnergyFunc, x) ::Real end

# for a smooth energy function V, we can calculate the work just
# by knowing "force" = grad(V), because dV/dt = grad(V) . dx/dt
function force(::SmoothEnergyFunc, F, x) ::Nothing end

# for a non-smooth V, grad(V) is not defined. Must implement
# the work directly
function work(::NonSmoothEnergyFunc, dxdt, x, sys, t) ::Real end

struct KuramotoEnergy <: SmoothEnergyFunc end;
struct MaxEnergy <: NonSmoothEnergyFunc end;

# energy function appropriate for systems where the x_i represent angles
# and we're concerned about phase-synchronization. x_i = x_j for all
# i, j is then the zero-energy state.
function energy(::KuramotoEnergy, x)
    n = length(x)
    sin_x = @~ sin.(x)
    cos_x = @~ cos.(x)
    # normalized between 0 and 1
    return 1.0 - (sum(sin_x)^2 + sum(cos_x)^2)/n^2
end

# force = gradient of energy w.r.t. the system states x
# force should be implemented as an in-place function, writing the output
# to a pre-allocated vector F
function force(::KuramotoEnergy, F, x)
    n = length(x)

    # do this in a crafty way (with lazyarrays) to avoid allocations
    sin_x = @~ sin.(x)
    cos_x = @~ cos.(x)
    sum_sin, sum_cos = sum(sin_x), sum(cos_x)
    @. F = -2.0 * (sum_sin*cos_x - sum_cos*sin_x)/n^2
end

@inline arc_bounds(x) = extrema(@~ mod2pi.(x))

@inline function arc_bound_indices(x)
    lb, ub = arc_bounds(x)
    all_idxs = 1:length(x)
    𝓘_max = Iterators.filter(i -> mod2pi(x[i]) == ub, all_idxs)
    𝓘_min = Iterators.filter(i -> mod2pi(x[i]) == lb, all_idxs)
    return 𝓘_min, 𝓘_max
end

geodesic(xi, xj) = mod(abs(xi-xj), pi)

function worst_geodesic(x)
    n = length(x)
    return maximum(geodesic(x[i], x[j]) for i=1:n for j=i+1:n)
end

energy(::MaxEnergy, x) = worst_geodesic(x)

function worst_case_pairs(x)
    V = worst_geodesic(x)
    n = length(x)
    return ((i, j) for i=1:n for j=i+1:n if geodesic(x[i], x[j]) == V)
end

function work(::MaxEnergy, dxdt::AbstractVector{<:Real}, x::AbstractVector{<:Real}, sys, t=0.0)
    rhs(dxdt, x, sys, t)
    return maximum(sign(x[i]-x[j])*(dxdt[i] - dxdt[j]) for (i, j) in worst_case_pairs(x))
    # dxdt_max = maximum(dxdt[i] for i in 𝓘_max)
    # dxdt_min = minimum(dxdt[i] for i in 𝓘_min)
    # return dxdt_max - dxdt_min
end

# work = total time derivative of energy (Lie derivative of energy). This works
# out to force ⋅ dxdt, so can automatically calculate for a generic EnergyFunc
# provided force has been defined. Note: this is in-place, requiring
# a pre-allocated vectors dxdt, F to store the force and rhs of the system.
function work(V::SmoothEnergyFunc, F::AbstractVector{<:Real}, dxdt::AbstractVector{<:Real}, 
              x::AbstractVector{<:Real}, sys, t=0.0; compute_F=true)
    if compute_F
        force(V, F, x)
    end
    rhs(dxdt, x, sys, t)
    return F ⋅ dxdt
end

# generic out-of-place version of force, for convenience
function force(V::SmoothEnergyFunc, x)
    F = similar(x)
    return force(V, F, x)
end

# generic out-of-place versions of the above, for convenience
function work(V::SmoothEnergyFunc, x::AbstractVector{<:Real}, sys, t=0.0)
    dxdt = similar(x)
    F = similar(x)
    return work(V, F, dxdt, x, sys, t)
end

function work(V::NonSmoothEnergyFunc, x::AbstractVector{<:Real}, sys, t=0.0)
    dxdt = similar(x)
    return work(V, dxdt, x, sys, t)
end