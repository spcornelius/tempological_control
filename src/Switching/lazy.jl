struct LazyWorkBasedSwitchingStrategy{T1<:EnergyFunc, T2<:WorkCache} <: WorkBasedSwitchingStrategy{T1, T2}
    V::T1
    cache::T2
end

function LazyWorkBasedSwitchingStrategy(V::EnergyFunc,
                                        dims::Union{Tuple, Int},
                                        m::Int)
     LazyWorkBasedSwitchingStrategy(V, WorkCache(dims, m))
end

function switch_condition(strategy::LazyWorkBasedSwitchingStrategy, ss, x, t)
    # when work in currently-active subsystem crosses zero, switch
    current_work(strategy, ss, x, t) >= 0
end
