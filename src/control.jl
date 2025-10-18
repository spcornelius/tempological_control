struct SwitchingSignal{T <: Real}
    t::Vector{T}
    i::Vector{Int}
    t_switch::Vector{T}
    transitions::Vector{Pair{Int, Int}}

    function SwitchingSignal(t_save::Vector{T}, i_save::Vector{Int}) where {T <: Real}
        n = length(i_save)

        # time indices where subsystem is different than in previous timestep
        switch_idx = findall(j -> i_save[j] != i_save[j+1], 1:n-1)
        t_switch = t_save[switch_idx]
        transitions = [i_save[j] => i_save[j+1] for j in switch_idx]

        # active subsystem at t=0
        i₀ = i_save[1]

        # sanity checks
        if ~isempty(transitions)
            @assert i₀ != transitions[1].second
            @assert i₀ == transitions[1].first
        end

        i = Vector{Int}([i₀])
        t = zeros(T, 1)
        for (tₛ, transition) in zip(t_switch, transitions)
            push!(t, tₛ)
            push!(t, tₛ)
            push!(i, transition.first)
            push!(i, transition.second)
        end
        push!(t, t_save[end])
        push!(i, i_save[end])
        new{T}(t, i, t_switch, transitions)
    end
end

struct ControlResult{T1<:AbstractODESolution, T2<:SwitchingSignal}
    success::Bool
    sol::T1
    signal::T2
end

# function barrier for type stability
function _solve(prob, integrator, cbs)
    return solve(prob, integrator; callback=cbs)
end

function control(ss, strategy, x₀, V;
                 t_max::TmaxType = 1000.0, integrator=Tsit5(), jac=nothing,
                 V_tol = 1e-8, Δt = 0.01, terminate_at_fixed_pt::Bool = true,
                 abstol = 1e-6, reltol = 1e-6,
                 cb = nothing) where {TmaxType}

    state_size(ss) == size(x₀) ||
        error("Dimensions of initial condition x₀ don't match " *
              "those of switched system")

    converged(x) = energy(V, x) < V_tol

    # main callback for switching
    switch_cb = PeriodicCallback(strategy, Δt, 
                                 save_positions=(true, true))

    # end integration early if energy falls below energy_tol
    convergence_cb = DiscreteCallback((x, t, integrator) -> converged(x),
                                      integrator -> terminate!(integrator))
    
    terminate_cb = terminate_at_fixed_pt ? TerminateSteadyState(abstol, reltol) : nothing

    t_save = Vector{TmaxType}()
    i_save = Vector{Int}()

    affect! = integrator -> begin
        push!(t_save, integrator.t)
        push!(i_save, integrator.p.i)
    end

    saving_cb =  DiscreteCallback((x, t, integrator) -> true, affect!;
                                  save_positions=(false, false))

    cbs = CallbackSet(switch_cb, convergence_cb, saving_cb, terminate_cb, cb)

    t_span = (zero(TmaxType), t_max)

    ode_func = ODEFunction{true}(rhs, jac=jac)
    prob = ODEProblem(ode_func, x₀, t_span, ss)
    sol = _solve(prob, integrator, cbs)

    success = converged(sol.u[end])
    return ControlResult(success, sol, SwitchingSignal(t_save, i_save))
end
