#__precompile__()

module TemporalSync

using OrdinaryDiffEq
using DiffEqBase: DECallback, AbstractODESolution
using DiffEqCallbacks
using LazyArrays
using LinearAlgebra
using Graphs
using SimpleWeightedGraphs
using SparseArrays
using SparseArrays: AbstractSparseMatrix

import Base: broadcastable

include("./config.jl")
include("./odesys.jl")
include("./kuramoto.jl")
include("./energy.jl")
include("./generate.jl")
include("./Switching/Switching.jl")
include("./control.jl")

export  # config.jl
        ProjectConfig,

        # odesys.jl
        ODESys,
        rhs,
        jac,
        state_size,

        # switchedsys.jl
        SwitchedSys,
        num_subsystems,

        # kuramoto.jl
        KuramotoSys,

        # generate.jl
        gen_kuramoto,

        # energy.jl
        energy,
        force,
        work,
        EnergyFunc,
        SmoothEnergyFunc,
        NonSmoothEnergyFunc,
        KuramotoEnergy,
        MaxEnergy,
        arc_bounds,
        arc_bound_indices,
        geodesic,
        worst_geodesic,
        worst_case_pairs,
        broadcastable,

        # control.jl
        control, SwitchingSignal, ControlResult
end
