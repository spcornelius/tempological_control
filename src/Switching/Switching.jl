include("./switchedsys.jl")
include("./strategy.jl")
include("./workbased.jl")
include("./lazy.jl")
include("./eager.jl")

export  # switchedsys.jl
        SwitchedSys,

        # strategy.jl
        SwitchingStrategy, switch_condition, choose_next_subsystem,
        initial_subsystem, current_work,

        # workbased.jl
        WorkBasedSwitchingStrategy, WorkCache, update_subsystem_works!,
        minimal_work_subsystem,

        # lazy.jl
        LazyWorkBasedSwitchingStrategy,

        # eager.jl
        EagerWorkBasedSwitchingStrategy
