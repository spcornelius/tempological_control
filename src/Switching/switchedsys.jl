mutable struct SwitchedSys{T} <: ODESys
    const subsystems::Vector{T}
    i::Int

    function SwitchedSys(subsystems::AbstractVector{T}, i::Int) where {T <: ODESys}
        N = length(subsystems)

        # make sure subsystem array isn't empty
        N > 0 ||
            error("SwitchedSystem instances must have at least one subsystem.")

        # make sure all subsystems have identical phase-space dimension        
        mapfoldl(state_size, ==, subsystems) ||
            error("All subsystems must have same phase space dimension.")

        new{T}(subsystems, i)
    end
end

SwitchedSys(subsystems::NTuple{<:Any, T}, i) where {T <: ODESys} =
    SwitchedSys(collect(subsystems), i)
    
SwitchedSys(subsystems) = SwitchedSys(subsystems, 1)

# create setter for SwitchedSys.i
set(ss::SwitchedSys, ::Val{T}, value) where {T} = setfield!(ss, T, value)

@inline function set(ss::SwitchedSys, ::Val{:i}, value::Integer)
    N = length(ss.subsystems)
    1 ≤ value ≤ N || error("Subsystem index i must be between 1 and $N.")
    setfield!(ss, :i, value)
end

Base.setproperty!(ss::SwitchedSys, field::Symbol, value) = set(ss, Val(field), value)

state_size(ss::SwitchedSys) = state_size(ss.subsystems[1])
num_subsystems(ss::SwitchedSys) = length(ss.subsystems)

function rhs(dxdt, x, ss::SwitchedSys, t=0.0)
    @inbounds rhs(dxdt, x, ss.subsystems[ss.i], t)
end

function jac(J, x, ss::SwitchedSys, t=0.0)
    @inbounds jac(J, x, ss.subsystems[ss.i], t)
end
