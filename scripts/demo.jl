using CairoMakie
using ColorSchemes
using Distributions
using FilePathsBase: /, PosixPath
using IterTools
using LaTeXStrings
using Makie
using Printf
using ProgressMeter
using Random
using Statistics
using TempologicalControl
using TempologicalControl: ProjectConfig
using ThreadsX

n = 10
q = 0.44
m_min = 1
m_max = 100
num_m = 20
num_trials = 500
V_tol = 1e-6
Δt = 0.05
t_max = 300.0
ps = collect(0.05:0.05:0.5)

k_avg = q * (n-1)

ms = unique(Int.(round.(
    10.0.^(range(log10(m_min), log10(m_max), num_m))
)))

set_theme!(ProjectConfig.base_theme)

colormap = cgrad(:algae, length(ps), categorical=true)

# get data
iter = IterTools.product(ps, ms, 1:num_trials)

progress = Progress(length(iter))

success = ThreadsX.map(iter) do (p, m, _)
    snapshots = [gen_kuramoto(n, k_avg, p) for _ in 1:m]
    sys = SwitchedSys(snapshots)
    V = KuramotoEnergy()
    strategy = EagerWorkBasedSwitchingStrategy(V, n, m)
    x₀ = rand(Uniform(-pi, pi), n)
    result = control(sys, strategy, x₀, V;
                     t_max = t_max, Δt = Δt, V_tol = V_tol)
    next!(progress)
    result.success
end

# success is a 3D array of dimension length(ps) × length(ms) × num_trials; average by m
success = dropdims(mean(success; dims=3); dims=3)

fig = Figure()
ax = Axis(fig[1,1], xscale=log10, xlabel=L"m", ylabel="success rate")

for (p, s, c) in zip(ps, eachrow(success), colormap)
    lines!(ax, ms, s, color=c, label=@sprintf("%.2f", p))
end

ax.xticks = [1, 10, 100]

axislegend(ax, L"p", position=:rb)

fig