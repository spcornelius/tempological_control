using ArgParse
using DataFrames
using LaTeXStrings
using FileIO
using FilePaths
using FilePathsBase: /
using IterTools
using KernelDensity
using PyCall
using PyPlot
using JLD2
using Statistics
using StatsBase
using SplitApplyCombine
using TemporalSync

include("./options.jl")
include("./plot_options.jl")

art3d = pyimport("mpl_toolkits.mplot3d.art3d")

parse_settings = ArgParseSettings()

@add_arg_table parse_settings begin
    "data_file"
        help = "path to .jld2 file containing simulation data"
        required = true
        arg_type = String
    "opts_file"
        help = "path to .toml file containing plot configuration"
        required = true
        arg_type = String
end

function fill_between_3d(ax, x, y, z; kwargs...)
    verts = collect(zip(x, y, z))
    append!(verts, [(x[end], y[end], 0.0)])
    append!(verts, [(x[1], y[1], 0.0)])
    return ax.add_collection3d(art3d.Poly3DCollection([verts]; kwargs...))
end

function success_vs_dist_fig(data_file, opts_file)
    opts = from_toml(PlotOptions, opts_file)

    gridspec = pyimport("matplotlib.gridspec")
    cm = pyimport("matplotlib.cm")
    cmap = cm.get_cmap(opts.cmap)
    slice(i,j) = pycall(pybuiltin("slice"), PyObject, i,j)
    rc("text",usetex="False")
    rc("font", size=opts.font_size)
    rc("xtick", labelsize=opts.tick_label_size) 
    rc("ytick", labelsize=opts.tick_label_size)
    rc("xtick.major", size=opts.major_tick_size)
    rc("xtick.minor", size=opts.minor_tick_size)
    rc("ytick.major", size=opts.major_tick_size)
    rc("ytick.minor", size=opts.minor_tick_size)
    rc("axes", labelsize=opts.axis_label_size, labelpad=opts.axis_label_pad)

    data = jldopen(string(data_file), "r")

    V = KuramotoEnergy()
    V₀_min = data["opts"].V₀_min
    V₀_max = data["opts"].V₀_max
    V₀_log_bins = range(log10(V₀_min), 
                        log10(V₀_max),
                        length=opts.V₀_bins+1)
    V₀_bins = 10.0.^V₀_log_bins

    # bin centers (in log space)
    V₀ = 10.0.^((V₀_log_bins[2:end] + V₀_log_bins[1:end-1])/2)

    m = opts.m
    n = opts.n

    X₀ = data["X₀/$n"]

    # get index of which V₀_bin each initial condition is in
    idx = [searchsortedfirst(V₀_bins, v) - 1 for v in energy.(V, X₀)]

    fig = figure(1, figsize=(opts.fig_width, opts.fig_height))
    plt.clf()
    gs = gridspec.GridSpec(nrows=1, ncols=3, figure=fig, width_ratios=opts.width_ratios)

    ax1 = fig.add_subplot(get(gs, 0))
    ax2 = fig.add_subplot(get(gs, 1))
    ax3 = fig.add_subplot(get(gs, 2), projection="3d", proj_type="persp")

    p_range = opts.p_range
    p_min = minimum(p_range)
    p_max = maximum(p_range)

    # each value of p gives one probability density
    # keep track of the maximum "z" value encountered, to properly scale z-axis
    z_max = 0.0

    x = collect(0.0:0.005:1.0)
    bins = range(0.0-1e-3, 1.0+1e-3, length=opts.panelB.bins)
    x_hist = (bins[2:end] .+ bins[1:end-1])./2

    # ax1.hist(energy.(V, X₀), bins=V₀_bins, color="lightgray", alpha=0.5, weights=ones(length(X₀))./length(X₀))
    for p in p_range
        s = data["success/$n/$m/$p"]

        # success rate by IC (averaged over all systems)
        s_by_θ₀ = dropdims(mean(s, dims=2), dims=2)

        # now average within each bin of V₀
        s_by_V₀ = combine(DataFrames.groupby(DataFrame(s=s_by_θ₀, idx=idx), :idx), :s => mean)[:,2]

        # success rate by system (averaged over all ICs)
        s_by_system = dropdims(mean(s, dims=1), dims=1)

        color = cmap((p - p_min)/(p_max - p_min))
        ax1.semilogx(V₀, s_by_V₀, color=color, marker=opts.panelA.marker, ms=opts.panelA.marker_size)

        ax2.semilogx([], [], color=color, marker=opts.legend.marker, ms=opts.legend.marker_size, 
                     linestyle="none", label="$p")

        ik = InterpKDE(kde(s_by_system, bandwidth=opts.panelB.bandwidth))
        z_max = max(z_max, maximum(pdf(ik, x)))

        if opts.panelB.mode == "density"
            z = pdf(ik, x)
            fill_between_3d(ax3, x, -p*ones(length(x)), z, color=color, ec="#262626",alpha=opts.panelB.alpha,
                            clip_on=false)
        elseif opts.panelB.mode == "hist"
            h = fit(Histogram, s_by_system, bins)
            z_max = max(z_max, h.weights)
            ax3.bar(x_hist, h.weights, width=1.0/opts.panelB.bins, zs=-p, zdir="y", color=color, 
                    ec=color, alpha=opts.panelB.alpha)
        else
            error("""Plot mode must be either "density" or "hist"; got "$opts.mode".""")
        end
    end

    for (axis, which) in IterTools.product(["x", "y"], ["major", "minor"])
        ax1.tick_params(axis=axis, which=which, pad=opts.panelA.tick_pad)
    end

    ax1.set_ylabel("success rate")
    ax1.set_xlabel(L"V_0", labelpad=opts.panelA.x_axis_label_pad)
    ax1.set_xlim(V₀_min, V₀_max)
    ax1.set_xticks(opts.panelA.xticks)

    ax3.view_init(elev=opts.panelB.elev, azim=opts.panelB.azim)
    ax3.set_ylim(-p_max, -p_min)
    ax3.set_zlim(0.0, z_max)
    ax3.set_yticks([])
    ax3.set_xlabel("success rate", labelpad=opts.panelB.x_axis_label_pad)
    ax3.set_zlabel("prob. density", labelpad=opts.panelB.z_axis_label_pad)
    for (axis, which) in IterTools.product(["x", "y", "z"], ["major", "minor"])
        ax3.tick_params(axis=axis, which=which, pad=opts.panelB.tick_pad)
    end
    ax3.dist = opts.panelB.dist

    # panel labels; use fig.text instead of ax.text to ensure alignment between panels
    fig.text(opts.panelA.panel_label_x, opts.panel_label_y, "(a)", fontsize=opts.panel_label_size)
    fig.text(opts.panelB.panel_label_x, opts.panel_label_y, "(b)", fontsize=opts.panel_label_size)

    ax2.legend(loc=(opts.legend.loc_x, opts.legend.loc_y), frameon=true, 
               labelspacing=opts.panelB.legend_label_spacing,
               handletextpad=opts.panelB.handle_text_pad)
    ax2.set_title(L"p", x=opts.legend.title_x, y=opts.legend.title_y, fontsize=opts.legend.title_size)
    ax2.axis("off")

    out_file = string(ProjectConfig.FIG_ROOT / opts.fig_file_name)
    plt.savefig(out_file, bbox_inches="tight")
    run(`pdfcrop $out_file $out_file`);
end

args = parse_args(ARGS, parse_settings)
success_vs_dist_fig(args["data_file"], args["opts_file"])
