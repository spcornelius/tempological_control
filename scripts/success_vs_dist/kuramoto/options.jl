using Configurations
    
@option "sim" struct SimOptions
    num_sys::Int64
    num_pts::Int64

    V₀_min::Float64
    V₀_max::Float64

    n_range::Vector{Int64}
    m_range::Vector{Int64}
    k_range::Vector{Float64}
    p_range::Vector{Float64}

    V_tol::Float64
    Δt::Float64
    t_max::Float64   
end

@option "sys" struct SysOptions
    slurm_mem_per_cpu::String
    # number of cores to use if not on slurm; 0 = use all but one
    local_cpu_cores::Int64
end

@option struct Options
    output_file_name::String
    sim::SimOptions
    sys::SysOptions
end