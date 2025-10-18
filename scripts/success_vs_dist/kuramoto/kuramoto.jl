using Distributed
using ClusterManagers
using Configurations
using Defer
using FileIO
using FilePaths: PosixPath
using FilePathsBase: /
using IterTools
using JLD2
using Printf
using ProgressMeter
using ParallelProgressMeter
using SPCUtil
using TemporalSync: ProjectConfig

# load configuration
const SRC_PATH = PosixPath(dirname(Base.source_path()))
include("options.jl")
opts = from_toml(Options, string(SRC_PATH / "config.toml"))

# utility function to determine whether an IO stream is a TTY or not
# (don't want to output progress bars if stderr is being recorded to a file, e.g. on SLURM)
is_logging(io) = isa(io, Base.TTY) == false || (get(ENV, "CI", nothing) == "true")

setup_workers()

# new scope and the following line are to ensure workers get cleaned
# up even if there is an exception, etc.
@scope begin
    @defer rmprocs(procs)

    @everywhere begin
        @eval using ArbNumerics, Distributions, Random, ProgressMeter, ParallelProgressMeter, TemporalSync        
        
        include("./bessel_ratio.jl")

        const DIST = Uniform(log10($(opts.sim.V₀_min)), log10($(opts.sim.V₀_max)))
        const V = KuramotoEnergy()

        function gen_x₀(n)
            # the expected circular mean (i.e., magnitude of the mean phasor) of a set of
            # Von Mises random variates is I₁(κ)/I₀(κ), so find the κ that corresponds
            # (in expectation) to the desired initial Kuramoto energy, V₀

            # need to use ArbFloats so that inverse_bessel_ratio is numerically stable
            # at the extremes; convert to normal Float64 later
            V₀_desired = ArbFloat(10.0^rand(DIST))
            κ = Float64(inverse_bessel_ratio(sqrt(1-V₀_desired)))
            while true
                x₀ = rand(VonMises(0, κ), n)

                # # ensure energy is within desired range
                # if $(opts.sim.V₀_min) <= energy(V, x₀) <= $(opts.sim.V₀_max)
                #     return x₀
                # end
                V₀ = energy(V, x₀)
                close_enough = abs((V₀ - V₀_desired)/V₀_desired) <= 0.05
                in_range = $(opts.sim.V₀_min) <= V₀ <= $(opts.sim.V₀_max)
                if close_enough && in_range
                    return x₀
                end
            end
        end

        function trial(x₀, sys)
            n = length(x₀)
            m = num_subsystems(sys)
            strategy = EagerWorkBasedSwitchingStrategy(V, n, m)
            control(sys, strategy, x₀, V;
                    t_max = $(opts.sim.t_max),
                    Δt = $(opts.sim.Δt),
                    V_tol = $(opts.sim.V_tol), full_output=false)
        end
    end

    jldopen(OUTPUT_FILE, "w") do data
        data["opts"] = opts.sim

        param_combs = IterTools.product(opts.sim.n_range, opts.sim.m_range, 
                                        opts.sim.k_range, opts.sim.p_range)

        # number of simulations per parameter
        num_sims = opts.sim.num_sys * opts.sim.num_pts

        # label each progress bar with parameters
        desc(n, m, k, p) = @sprintf "n=%i; m=%i; p=%.2f: " n m p

        pbars =  map(p -> Progress(num_sims, desc=desc(p...)), param_combs)

        # need to flatten in order to provide to MultipleProgress
        pbars = reduce(vcat, pbars)

        progress = MultipleProgress(pbars, Progress(num_sims*length(param_combs), desc="All :"),
                                    enabled=!is_logging(stderr))

        println("Generating initial conditions...")
        for n in opts.sim.n_range
            X₀ = [gen_x₀(n) for _ in 1:opts.sim.num_pts]
            data["X₀/$n"] = X₀
        end

        println("Starting simulations...")
        @sync begin
            lk = ReentrantLock()
            for (i, params) in enumerate(param_combs)
                @async begin
                    n, m, k, p = params
                    systems = [SwitchedSys([gen_kuramoto(n, k, p) for _ in 1:m])
                               for _ in 1:opts.sim.num_sys]

                    # save systems
                    # JLD2 not thread-safe; make sure only one task writes at a time
                    lock(lk) do
                        data["systems/$n/$m/$p"] = systems
                    end    

                    X₀ = data["X₀/$n"]
                    success = pmap(IterTools.product(X₀, systems)) do (x₀,sys)
                        s = trial(x₀,sys)
                        next!(progress[i])
                        s
                    end
                    
                    # save results
                    lock(lk) do
                        data["success/$n/$m/$p"] = success
                    end
                end 
            end 
        end
    end
end

# move output file from scratch to final location
if IN_SLURM
    mv(OUTPUT_FILE, OUTPUT_FILE_DEST, force=true)
end

if @isdefined tmp_dir
    rm(tmp_dir, recursive=true, force=true)
end