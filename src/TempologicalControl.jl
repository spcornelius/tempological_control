#__precompile__()

module TempologicalControl

using OrdinaryDiffEq
using DiffEqBase
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
        KuramotoEnergy,
        broadcastable,

        # control.jl
        control, SwitchingSignal, ControlResult
end
