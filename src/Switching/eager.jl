struct EagerWorkBasedSwitchingStrategy{T1<:EnergyFunc, T2<:WorkCache} <: WorkBasedSwitchingStrategy{T1, T2}
    V::T1
    cache::T2
end

function EagerWorkBasedSwitchingStrategy(V::EnergyFunc,
                                         dims::Union{Tuple, Int},
                                         m::Int)
    EagerWorkBasedSwitchingStrategy(V, WorkCache(dims, m))
end

function switch_condition(strategy::EagerWorkBasedSwitchingStrategy, ss, x, t)
    # when there is a better subsystem, switch
    j = minimal_work_subsystem(strategy, ss, x, t; exclude_current=true)
    i = ss.i
    W = strategy.cache.W
    return W[j] < W[i]
end
