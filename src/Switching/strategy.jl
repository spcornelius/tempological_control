abstract type SwitchingStrategy end;

function switch_condition(strategy::SwitchingStrategy, ss, x, t) ::Bool end
function choose_subsystem(strategy::SwitchingStrategy, ss, x, t) ::Int end;
function initial_subsystem(strategy::SwitchingStrategy, ss, x, t) ::Int end;

function DiffEqBase.DiscreteCallback(strategy::SwitchingStrategy; kw...)
    condition = (x, t, integrator) -> switch_condition(strategy, integrator.p, x, t)
    affect! = integrator -> integrator.p.i = choose_subsystem(strategy, integrator.p, integrator.u, integrator.t)
    initialize = (cb, x, t, integrator) -> integrator.p.i = initial_subsystem(strategy, integrator.p, integrator.u, integrator.t)

    return DiscreteCallback(condition, affect!;
                            initialize=initialize,
                            kw...)
end

function DiffEqCallbacks.PeriodicCallback(strategy::SwitchingStrategy, Δt::Number; kw...)
    affect! = integrator -> integrator.p.i = choose_subsystem(strategy, integrator.p, integrator.u, integrator.t)
    return PeriodicCallback(affect!, Δt, initial_affect=true; kw...)
end
