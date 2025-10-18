using CairoMakie
using Colors
using Distributions
using DrWatson
using EmpiricalCDFs
using Fontconfig
using Formatting
using FreeTypeAbstraction
using FromFile
using Graphs
using Infiltrator
using LinearAlgebra
using Makie
using MathTeXEngine
using Random
using ProgressMeter
using SparseArrays
using StatsBase
using TemporalSync
using TemporalSync: ProjectConfig
using ThreadsX
using UnPack

DrWatson.@quickactivate

@from "./util.jl" import Directed, Undirected, work_mean_and_std, all_works
@from "./gumbel_util.jl" import ReverseGumbel, gumbel_loc_and_scale, QuantileExact, 
                                QuantileAsymptotic, WikipediaExact, WikipediaAsymptotic

const SUB_DIR = "theory"
const DATA_PATH = joinpath(datadir(), SUB_DIR)

mode = Undirected()

n = 100
k_avg = 6.0
ρs_hists = [0.6, 0.7, 0.8, 0.9]
ms_hists = [5, 500]
ρs_lines = [0.6, 0.7, 0.8, 0.9]
ms_lines = unique(Int.(ceil.(logrange(2, 20000, length=20))))
trials = 1_000
dθ = 0.05

# ranges for bins of all histograms (want to ensure bins are equal width across histograms)
w_min = -0.1
w_max = 0.1
bins = 50
bin_edges = range(w_min, w_max, bins + 1)
bin_centers = (bin_edges[2:end] .+ bin_edges[1:end-1])./2
neg_idxs = bin_centers .< 0
nonneg_idxs = bin_centers .>= 0

x0 = zeros(n)
for i in 2:n
    x0[i] = x0[i-1] + dθ
end

@show energy(KuramotoEnergy(), x0)

function lighten(c, factor)
    c = HSV(c)
    @unpack h, s, v = c
    s = s * factor
    return HSV{Float32}(h, s, v)
end

#colors = Dict(ρ=>Makie.wong_colors()[i] for (i, ρ) in enumerate(ρs_lines))
saturated_colors = [colorant"rgb(04,146,34)", colorant"rgb(31, 178, 196)", colorant"rgb(27, 98, 165)", colorant"rgb(129,78,175)"]
saturated_colors = Dict(ρ=>c for (ρ, c) in zip(ρs_hists, saturated_colors))

unsaturated_colors = [colorant"rgb(136,220,120)", colorant"rgb(142,211,222)", colorant"rgb(158,186,226)", colorant"rgb(184,158,204)"]
unsaturated_colors = [lighten(c, 0.85) for c in unsaturated_colors]
unsaturated_colors = Dict(ρ=>c for (ρ, c) in zip(ρs_hists, unsaturated_colors))

colors = Dict(ρ=>c for (ρ, c) in zip(ρs_lines, saturated_colors))
function font_file(pattern::Fontconfig.Pattern)
    Fontconfig.format(Fontconfig.match(pattern), "%{file}")
end

font_file(str::String) = font_file(Fontconfig.Pattern(str))

d = set_texfont_family!(regular = font_file("arial"))
set_texfont_family!()
MathTeXEngine.default_font_families["Arial"] = d

set_theme!(ProjectConfig.base_theme)

update_theme!(
    Theme(
        #size = (450, 175),
        size = (333, 133),
        fontsize = 8,
        figure_padding = (4, 10, 0, 2),
        Axis = (
            titlesize = 10,
            xlabelsize = 10,
            ylabelsize = 10,
            xlabelpadding = 3,
            ylabelpadding = 8,
            topspinevisible=false, 
            rightspinevisible=false,
            spinewidth=0.75,
        ),
        Label = (
            #font = "Arial",
            fontsize = 12,
            #color = ALMOST_BLACK,
            #halign = :left
        )
    )
)

#x0[2:2:end] .= pi/2
#x0 = rand(Uniform(-pi, pi), n)
#x0 = pi*collect(0:n-1)./n

function rand_weights(k, ρ)
    w = rand(Uniform(), k)
    for i in eachindex(w)
        if rand() < ρ
            w[i] *= -1
        end
    end
    return w
end

function gen_fast(mode::Directed, n, k_avg, ρ; kw...)
    p = k_avg/n
    A = sprand(n, n, p, L -> rand_weights(L, ρ))
    for i in 1:n
        A[i, i] = 0.0
    end
    return KuramotoSys(A)
end

function gen_fast(mode::Undirected, n, k_avg, ρ; kw...)
    p = k_avg/n/2
    A = sprand(n, n, p, L -> rand_weights(L, ρ))
    
    # make symmetric
    A = triu(A)
    A = A + A'

    for i in 1:n
        A[i, i] = 0.0
    end
    return KuramotoSys(A)
end

function work_trial(mode, n, k_avg, ρ, m, x0)
    systems = (gen_fast(mode, n, k_avg, ρ; force_connected=false, force_unstable=false) for _ in 1:m)
    return minimum(all_works(systems, x0))
end

function work_sim(config)
    @unpack mode, n, k_avg, ρ, m, x0, trials = config
    p = Progress(trials)
    best_works = ThreadsX.mapi(1:trials) do _
        result = work_trial(mode, n, k_avg, ρ, m, x0)
        next!(p)
        result
    end
    return @strdict(best_works)
end

function success_trial(mode, n, k_avg, ρ, m, x0)
    # this can be much faster than the above as the "any" will bail
    # as soon as there is a success. Downside is we don't get the minimal actual work...
    systems = (gen_fast(mode, n, k_avg, ρ; force_connected=false, force_unstable=false) for _ in 1:m)
    V = KuramotoEnergy()
    F = similar(x0)
    dxdt = similar(x0)
    force(V, F, x0)
    return any(s -> work(V, F, dxdt, x0, s; compute_F=false) < 0, systems)
end

fig = Figure()

#ylabel1 = "density"
ax1 = Axis(fig[1,1], title="$(Formatting.format(ms_hists[1],commas=true)) snapshots")

xlabel2 = L"{\fontfamily{Arial}minimal work,} $W_\mathit{\min}$"
ax2 = Axis(fig[2,1], xlabel=xlabel2, title="$(Formatting.format(ms_hists[2],commas=true)) snapshots")

linkxaxes!(ax1, ax2)
linkyaxes!(ax1, ax2)

xlabel2 = L"{\fontfamily{Arial}number of snapshots,} $m$"
ylabel2 = L"\mathbb{P}(W_\min\;<\;0)"
ax3 = Axis(fig[:,2], xscale=log10, ylabel=ylabel2, xlabel=xlabel2)

lins = Any[]
for ρ in ρs_lines
    μ, σ = work_mean_and_std(mode, x0, n, k_avg, ρ)
    @show μ, σ, μ/σ

    pct_negative = Float64[]
    best_works = Vector{Float64}[]

    for m in ms_lines
        config = @dict(mode, n, k_avg, ρ, m, x0, trials)
        #config = Dict("n" => n, "k_avg" => k_avg, "ρ" => ρ, "m" => m, "x0" => x0, "trials" => trials)
        data, _ = produce_or_load(work_sim, config, DATA_PATH; filename=hash, verbose=false, tag=false)
        bw = data["best_works"]
        num_negative = count(w -> w < 0, bw)
        push!(best_works, bw)
        push!(pct_negative, num_negative/trials)
    end

    theory = Float64[]
    m_range = logrange(minimum(ms_lines), maximum(ms_lines), 100)
    for m in m_range
        b, a = gumbel_loc_and_scale(m, QuantileExact())
        loc = σ* b - μ
        scale = σ * a
        d = ReverseGumbel(loc, scale)
        #push!(theory, reverse_gumbel_cdf(0, a, b, μ, σ))
        push!(theory, cdf(d, 0.0))
    end

    c = saturated_colors[ρ]
    l = lines!(ax3, m_range, theory, label="$ρ", color=c, linewidth=1.5)
    push!(lins, l)
    scatter!(ax3, ms_lines, pct_negative, color=c, markersize=7)

    # hist!(ax1, best_works[4], bins=30, normalization=:pdf, color=colors[ρ])
    # hist!(ax1, best_works[15], bins=30, normalization=:pdf, color=colors[ρ])
end

x = range(w_min, w_max, 1000)
for ρ in ρs_hists
    μ, σ = work_mean_and_std(mode, x0, n, k_avg, ρ)
    for (ax, m) in zip((ax1, ax2), ms_hists)
        config = @dict(mode, n, k_avg, ρ, m, x0, trials)
        data, _ = produce_or_load(work_sim, config, DATA_PATH; filename=hash, verbose=false, tag=false)
        bw = data["best_works"]
        hist!(ax, bw, color=unsaturated_colors[ρ], bins=range(w_min, w_max, 40), normalization=:pdf)

        b, a = gumbel_loc_and_scale(m, QuantileExact())
        loc = σ* b - μ
        scale = σ * a
        d = ReverseGumbel(loc, scale)
        lines!(ax, x, pdf.(d, x), color=saturated_colors[ρ], linewidth=1.25)
        # @show maximum(pdf.(d, x))
        @show (ρ, cdf(d, 0.0))
    end
end

ylims!(ax2, 0.0, 50.0)
        # leg = axislegend(ax, L"p", position=:rb, labelsize=8, rowgap=0)
elems = [MarkerElement(marker=:rect, color=unsaturated_colors[ρ], strokecolor=saturated_colors[ρ], strokewidth=1) for ρ in ρs_lines]
# Legend(fig[:,2], lins, ["0.6" , "0.7", "0.8", "0.9"], L"p", tellheight=false, tellwidth=false, halign=:right, valign=:bottom,
#        labelsize=8, rowgap=-9, titlegap=0)
Legend(fig[:,2], elems, ["0.6" , "0.7", "0.8", "0.9"], L"\;\;\;\;\;p", tellheight=false, tellwidth=false, halign=:right, valign=:bottom,
       labelsize=6, rowgap=-9, titlegap=0, patchlabelgap=-2, padding=(0, 0, 0, 0))

xlims!(ax2, w_min, w_max)
hidexdecorations!(ax1)   
hideydecorations!(ax1)
hideydecorations!(ax2)
xlims!(ax3, 1, 1.3*maximum(ms_lines))

Label(fig[1, 1, TopLeft()], "(a)", padding=(0, 0, 4, 0))
Label(fig[2, 1, TopLeft()], "(b)", padding=(0, 0, 4, 0))
Label(fig[1, 2, TopLeft()], "(c)", padding=(0, 0, 4, 0))

rowgap!(fig.layout, 1, Relative(0.05))
colgap!(fig.layout, 1, Relative(0.04))

# Translate spines
x1, x2 = ax1.xaxis.attributes.endpoints[]
xl, xo = x1
xr = x2[1]
y1, y2 = ax1.yaxis.attributes.endpoints[]
yo, yb = y1
yt = y2[2]
new_xo = yb + (yt - yb)/2.0
new_yo = xl + (xr - xl)/2.0

ax1.xaxis.attributes.endpoints[] = ([xl, new_xo], [xr, new_xo])
ax1.yaxis.attributes.endpoints[] = ([new_yo, yb], [new_yo, yt])

x1, x2 = ax2.xaxis.attributes.endpoints[]
xl, xo = x1
xr = x2[1]
y1, y2 = ax2.yaxis.attributes.endpoints[]
yo, yb = y1
yt = y2[2]
new_xo = yb + (yt - yb)/2.0
new_yo = xl + (xr - xl)/2.0

ax2.xaxis.attributes.endpoints[] = ([xl, new_xo], [xr, new_xo])
ax2.yaxis.attributes.endpoints[] = ([new_yo, yb], [new_yo, yt])

ax1.yaxis.elements[:axisline].linewidth = 1.5
ax2.yaxis.elements[:axisline].linewidth = 1.5

#text!(ax1, -0.01, 40; text=L"$V${\fontfamily{Arial} decreasing}", align=(:right, :baseline))
#text!(ax1, 0.01, 40; text=L"$V${\fontfamily{Arial} increasing}", align=(:left, :baseline))

save_file = joinpath(plotsdir(), "theory.pdf")
save(save_file, fig, pt_per_unit=1)
display(fig)
