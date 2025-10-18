struct WorkCache{R<:Real, D}
    W::Vector{R}
    F::Array{R, D}
    dxdt::Array{R, D}

    function WorkCache(R::Type, dims::Union{Tuple,Int}, m::Int)
        W = zeros(R, m)
        F = zeros(R, dims)
        dxdt = zeros(R, dims)
        new{R, length(dims)}(W, F, dxdt)
    end
end

WorkCache(dims, m) = WorkCache(Float64, dims, m)

abstract type WorkBasedSwitchingStrategy{T1<:EnergyFunc, T2<:WorkCache} <: SwitchingStrategy end;

function current_work(strat::WorkBasedSwitchingStrategy, ss, x, t=0.0)
    cache = strat.cache
    return work(strat.V, cache.F, cache.dxdt, x, ss, t; compute_F=true)
end

function update_subsystem_works!(strat::WorkBasedSwitchingStrategy{<:SmoothEnergyFunc, <:Any},
                                 ss, x, t=0.0)
    cache = strat.cache
    m = num_subsystems(ss)
    force(strat.V, cache.F, x)
    for i=1:m
        cache.W[i] = work(strat.V, cache.F, cache.dxdt, x,
                          ss.subsystems[i], t; compute_F=false)
    end
end

function update_subsystem_works!(strat::WorkBasedSwitchingStrategy{<:NonSmoothEnergyFunc, <:Any},
                                 ss, x, t=0.0)
    cache = strat.cache
    m = num_subsystems(ss)
    for i=1:m
        cache.W[i] = work(strat.V, cache.dxdt, x,
                          ss.subsystems[i], t)
    end
end

function minimal_work_subsystem(strat::WorkBasedSwitchingStrategy,
                                ss, x, t=0.0;
                                exclude_current=false)
    update_subsystem_works!(strat, ss, x, t)
    W = strat.cache.W
    tmp = W[ss.i]
    if exclude_current
        W[ss.i] = Inf
    end
    i_min = argmin(W)
    W[ss.i] = tmp
    return i_min
end

function choose_subsystem(strategy::WorkBasedSwitchingStrategy, ss, x, t)
    minimal_work_subsystem(strategy, ss, x, t; exclude_current=false)
end

function initial_subsystem(strategy::WorkBasedSwitchingStrategy, ss, x, t)
    minimal_work_subsystem(strategy, ss, x, t; exclude_current=false)
end
