using CairoMakie
using Colors
using Configurations
using FilePaths: PosixPath
using FilePathsBase: /
using FromFile
using Graphs
using Infiltrator
using JLD2
using LaTeXStrings
using LinearAlgebra
using Makie
using SparseArrays
using SPCUtil
using TemporalSync
using YAML

CairoMakie.activate!()

this_dir = PosixPath(@__DIR__())

@from "schematic_options.jl" import SchematicOptions, schematic_theme
@from "schematic_util.jl" import net_schematic_fig
@from "vectorfield.jl" import vectorfield!

# Makie theme-ing
set_theme!(ProjectConfig.base_theme)
update_theme!(schematic_theme)

# load configuration
config_file_path = string(this_dir / "schematic_config.yml")
config_data = YAML.load_file(config_file_path; dicttype=Dict{String, Any})
opts = from_dict(SchematicOptions, config_data)

# common stuff
const ALMOST_BLACK = colorant"#262626"
const V = KuramotoEnergy()
const strategy = EagerWorkBasedSwitchingStrategy(V, 3, 2)

const tick_locs = [-pi, 0, pi]
const tick_labels = [L"-\pi", L"0", L"+\pi"]
const ticks = (tick_locs, tick_labels)
pp_axis(loc) = Axis(loc, xlabel=L"\theta_1", ylabel=L"\theta_2", 
                    xticks=ticks, yticks=ticks, aspect=1)

snapshot_colors = [colorant"#4daf4a",
                   colorant"rgb(250,184,10)"]#colorant"#984ea3"]

##############################################################################
# vector fields for snapshots
##############################################################################
function panel_a(grid, g1, g2, opts)
    net_schematic_fig(grid[1, 1], g1, opts)
    net_schematic_fig(grid[1, 2], g2, opts)

    pp1_ax = pp_axis(grid[2, 1])  
    pp2_ax = pp_axis(grid[2, 2])
    hideydecorations!(pp2_ax)
    pp1_ax.ylabelpadding = 4

    v1 = vectorfield!(pp1_ax, g1)
    v2 = vectorfield!(pp2_ax, g2)

    # set both color ranges to same scale, symmetric around W = 0
    # this is a bit kludge-y; the contourf plot is the 1st plot
    # within the 1st plot of the result of vectorfield!
    c1 = v1.plots[1].plots[1]
    c2 = v2.plots[1].plots[1]
    W_mag = maximum(abs.([c1.colorrange[]..., 
                          c2.colorrange[]...]))
    colorrange = (-W_mag, W_mag)
    c1.colorrange = colorrange
    c2.colorrange = colorrange

    for ax in [pp1_ax, pp2_ax]
        xlims!(ax, -pi, pi)
        ylims!(ax, -pi, pi)
        hidespines!(ax)
        ax.xgridvisible = false
        ax.ygridvisible = false
    end
    nothing
end

##############################################################################
# successful trajectory for switched system
##############################################################################
function panel_b(grid, u, t, signal, opts)
    # get θ₁ and θ₂ relative to θ₃
    θ₁ = @view u[1, :]
    θ₂ = @view u[2, :]

    # cumulative arc length
    s = cumsum(sqrt.(diff(θ₁).^2 + diff(θ₂).^2))

    # express as a fraction between 0 and 1
    s ./= s[end]
    pushfirst!(s, 0.0)

    # arrow locations (x, y)
    I_arrow = searchsortedfirst.(Ref(s), opts.arrow_locs)
    θ_arrow = u[1:2, I_arrow]
    #σ_arrow = signal.i[I_arrow[:]]

    # directions of arrows (u, v)
    # approximate as the 2nd-next state in the time series,
    # assuming it's been saved at sufficiently high time resolution
    θ_arrow_next = u[1:2, I_arrow .+ 2]

    arrow_dirs = (θ_arrow_next .- θ_arrow)'
    arrow_dirs = Point2f.(eachrow(arrow_dirs))
    arrow_locs = Point2f.(eachrow(θ_arrow'))

    pp_switched_ax = pp_axis(grid[1, 1])
    line_colors = snapshot_colors[signal.i[searchsortedfirst.(Ref(signal.t), t)]]
    arrow_colors = snapshot_colors[signal.i[searchsortedfirst.(Ref(signal.t), t[I_arrow])]]

    vlines!(pp_switched_ax, 0, color=:lightgray, linewidth=0.5)
    hlines!(pp_switched_ax, 0, color=:lightgray, linewidth=0.5)

    lines!(pp_switched_ax, θ₁, θ₂,  linewidth=1.75, color=line_colors)
    #lines!(pp_switched_ax, θ₁, θ₂,  linewidth=1, color=:black)
    scatter!(pp_switched_ax, [Point2f(u[1:2, 1])], color=:black, markersize=6)
    arrows!(pp_switched_ax, arrow_locs, arrow_dirs, arrowsize=18, color=arrow_colors)
    xlims!(pp_switched_ax, -pi, pi)
    ylims!(pp_switched_ax, -pi, pi)



    #hidespines!(pp_switched_ax)
    pp_switched_ax.xticksvisible = false
    pp_switched_ax.yticksvisible = false
    pp_switched_ax.xticklabelsvisible = false
    pp_switched_ax.yticklabelsvisible = false
end

##############################################################################
# energy, work, and switching signal over time
##############################################################################
function panel_c(grid, sys, t, u, signal, opts)
    energy_ax = Axis(grid[1, 1], ylabel=L"V", yscale=log10)
    work_ax = Axis(grid[2, 1], ylabel=L"\frac{W}{V}")
    ss_ax = Axis(grid[3, 1], ylabel=L"\sigma", xlabel=L"t")

    i_max = findfirst(u_ -> energy(V, u_) < opts.V_thresh, eachcol(u))
    t_max = t[i_max]
    for ax in (energy_ax, work_ax, ss_ax)
        ax.topspinevisible = false
        ax.rightspinevisible = false
        xlims!(ax, 0.0, t_max)
        #ax.yaxis.elements[:labeltext].attributes[:rotation] = 0.0
        #ax.yaxis.elements[:labeltext].attributes[:align] = (:center, :center)
    end

    ylims!(energy_ax, opts.V_thresh, nothing)

    energy_ax.yticks = LogTicks(WilkinsonTicks(3))
    ss_ax.xticks = WilkinsonTicks(5, k_min=5, k_max=5)

    V_ = energy.(V, eachcol(u))
    lines!(energy_ax, t, V_, linewidth=1, color=ALMOST_BLACK)

    W1 = work.(V, eachcol(u), sys.subsystems[1])
    W2 = work.(V, eachcol(u), sys.subsystems[2])
    # lines!(work_ax, t, W1./V_, color=colorant"#1b9e77", linewidth=1, label="1")
    # lines!(work_ax, t, W2./V_, color=colorant"#d95f02", linewidth=1, label="2")
    # lines!(work_ax, t, W1./V_, color=colorant"#E1AB2C", linewidth=1, label="1")
    # lines!(work_ax, t, W2./V_, color=colorant"#33AF7A", linewidth=1, label="2")
    hlines!(work_ax, 0, color=ALMOST_BLACK, linewidth=0.5)
    lines!(work_ax, t, W1./V_, color=snapshot_colors[1], linewidth=1.25, label="1")
    lines!(work_ax, t, W2./V_, color=snapshot_colors[2], linewidth=1.25, label="2")
    work_ax.bottomspinevisible=false
    #bf3c1a

    lines!(ss_ax, signal.t, signal.i, color=:black)

    ylims!(ss_ax, 0.9, 2.1)
    ss_ax.yticks = [1, 2]
    hidexdecorations!(work_ax)
    hidexdecorations!(energy_ax)

    energy_ax.yticklabelspace = tight_yticklabel_spacing!(energy_ax)
    align_ylabels!(energy_ax, work_ax, ss_ax)
    nothing
end

function schematic_fig(opts)
    fig = Figure(resolution=(opts.total_width, opts.total_height),
                 figure_padding=(-4, 1, 0, 20))

    ##############################################################################
    # panel A: network diagrams and phase planes
    ##############################################################################
    gridA = fig[1, 1] = GridLayout()

    # alignmode = Mixed(
    #     #left = Makie.Protrusion(0),
    #     #right = Makie.Protrusion(0),
    #     bottom = 0,
    #     #top = Makie.Protrusion(35),
    # )
    Box(gridA[1, 1], linestyle = :dash, 
        valign=:top, 
        height=67, 
        alignmode=Outside(0, 0, 0, -4),
        color=:white,
        cornerradius=7,
        strokewidth=2.0,
        strokecolor=snapshot_colors[1]
        )

    Box(gridA[1, 2], linestyle = :dash, 
        valign=:top, 
        height=67, 
        alignmode=Outside(0, 0, 0, -4),
        color=:white,
        cornerradius=7,
        strokewidth=2.0,
        strokecolor=snapshot_colors[2]
        )

    Label(gridA[1, 1], " Snapshot 1", color=snapshot_colors[1],
          alignmode=Outside(0, 0, 0, -65), tellwidth=false, tellheight=true)

    Label(gridA[1, 2], " Snapshot 2", color=snapshot_colors[2],
          alignmode=Outside(0, 0, 0, -65), tellwidth=false, tellheight=true)

    data_file = ProjectConfig.DATA_ROOT / opts.data_file
    data = jldopen(string(data_file), "r")

    g1, g2 = data["graph_pairs"][opts.index]
    θ₀ = data["θ₀"][opts.index]
    θ₀ = [-0.0, 2.0, 0.0]

    panel_a(gridA, g1, g2, opts.net)
    Label(gridA[1, 1, TopLeft()], "(a)", padding=(5, 0, -10, -30))

    ##############################################################################
    # panel B: switched trajectory
    ##############################################################################
    sys = SwitchedSys(KuramotoSys.([g1, g2]))
    res = control(sys, strategy, θ₀, V; Δt=0.01, t_max=100.0, V_tol=1e-3)
    sol = res.sol

    # measure θ₁ and θ₂ relative to θ₃
    t = sol.t
    u = sol[:,:]
    u = (u' .- u'[:, 3])'

    gridB = fig[1, 2] = GridLayout()
    panel_b(gridB, u, t, res.signal, opts)
    Label(gridB[1, 1, TopLeft()], "(b)", padding=(0, -10, -10, -30))

    ##############################################################################
    # panel C: energy, work, and switching signal over time
    ##############################################################################
    gridC = fig[1, 3] = GridLayout()
    panel_c(gridC, sys, t, u, res.signal, opts)
    Label(gridC[1, 1, TopLeft()], "(c)", padding=(15, 0, -10, -30))

    ##############################################################################
    # Clean-up and layout
    ##############################################################################

    # un-rotate all y-axis tick_labels
    for child in fig.content
        if isa(child, Axis)
            child.yaxis.elements[:labeltext].attributes[:rotation] = 0.0
            child.yaxis.elements[:labeltext].attributes[:align] = (:center, :center)
        end
    end

    for (i, gap) in enumerate(opts.outer_colgaps)
        colgap!(fig.layout, i, Relative(gap))
    end

    for (i, size) in enumerate(opts.outer_colsizes)
        colsize!(fig.layout, i, Relative(size))
    end

    rowgap!(gridA, Relative(opts.panel_a_rowgap))
    for (i, size) in enumerate(opts.panel_a_rowsizes)
        rowsize!(gridA, i, Relative(size))
    end

    rowgap!(gridC, Relative(opts.panel_c_rowgap))
    for (i, size) in enumerate(opts.panel_c_rowsizes)
        rowsize!(gridC, i, Relative(size))
    end

    resize_to_layout!(fig)
    Makie.trim!(fig.layout)

    return fig
end

for i in 53:53
    opts.index = i
    fig = schematic_fig(opts)
    save_file = ProjectConfig.FIG_ROOT / opts.save_file
    save(string(save_file), fig; pt_per_unit=0.75)
    display(fig)
end