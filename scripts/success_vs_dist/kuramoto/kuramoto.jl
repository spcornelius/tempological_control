using Distributions
using Configurations
using Defer
using FileIO
using FilePathsBase: /, PosixPath
using FromFile
using IterTools
using JLD2
using Printf
using ProgressMeter
using Random
using TempologicalControl
using TempologicalControl: ProjectConfig
using ThreadsX

# load configuration
@from "./options.jl" import Options
src_root = PosixPath(dirname(Base.source_path()))
opts = from_toml(Options, string(src_root / "config.toml"))

output_file = string(ProjectConfig.DATA_ROOT / "success_vs_dist" / "kuramoto" / opts.output_file_name)

function gen_x₀(n, opts)
    # distribution of desired initial energies
    V₀_dist = Uniform(log10(opts.V₀_min), log10(opts.V₀_max))
    V₀_desired = 10.0^rand(V₀_dist)

    while true
        x₀ = rand_ic_with_energy(n, V₀_desired; rtol=opts.rtol)
        V₀ = energy(KuramotoEnergy(), x₀)
        if opts.V₀_min <= V₀ <= opts.V₀_max
            return x₀
        end
    end
end

function trial(x₀, sys, opts)
    n = length(x₀)
    m = num_subsystems(sys)
    V = KuramotoEnergy()
    strategy = EagerWorkBasedSwitchingStrategy(V, n, m)
    control(sys, strategy, x₀, V;
            t_max = opts.t_max,
            Δt = opts.Δt,
            V_tol = opts.V_tol)
end

function main(output_file, opts)
    jldopen(output_file, "w") do data
        data["opts"] = opts

        param_combs = IterTools.product(opts.n_range, opts.m_range, 
                                        opts.k_range, opts.p_range)

        # number of simulations per parameter
        num_sims = opts.num_sys * opts.num_pts

        progress = Progress(num_sims*length(param_combs), desc="All :", enabled=!is_logging(stderr))

        println("Generating initial conditions...")
        for n in opts.n_range
            X₀ = [gen_x₀(n, opts) for _ in 1:opts.num_pts]
            data["X₀/$n"] = X₀
        end

        for (i, params) in enumerate(param_combs)
            n, m, k, p = params
            
            subsystem_sets = [[gen_kuramoto(n, k, p) for _ in 1:m] for _ in 1:opts.num_sys]
            snapshot_sets = [[s.A for s in subsystem_set] for subsystem_set in subsystem_sets] 

            # save MATRICES to avoid putting custom types into JLD2
            data["snapshots/$n/$m/$p"] = snapshot_sets
            X₀ = data["X₀/$n"]

            success = ThreadsX.map(IterTools.product(X₀, subsystem_sets)) do (x₀, subsystem_set)
                sys = SwitchedSystem(subsystem_set)
                s = trial(x₀, sys, opts)
                next!(progress)
                s.success
            end
            
            # save results
            data["success/$n/$m/$p"] = success
        end 
    end
end

if isfile(output_file)
    while true
        print("Output file $output_file exists. Overwrite (y/N)? ")
        response = lowercase(strip(chomp(readline())))
        if isempty(response) || response == "n"
            exit()
        elseif response == "y"
            break
        end
    end
end

main(output_file, opts.sim)
