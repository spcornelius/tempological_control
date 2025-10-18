using Colors
using Configurations
import Configurations: from_dict
using FromFile
using Makie

# custom type converter for matrices; should be specified in YAML as, for example:
# I = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
function Configurations.from_dict(::Type{OptionType}, of::OptionField, ::Type{T}, x::V) where 
        {OptionType, T<:AbstractMatrix, V<:AbstractVector{<:AbstractVector{<:Any}}}
    return convert(T, reduce(hcat, x)')
end

@option struct NetworkOptions
    node_color::String
    text_color::String
    pos_link_color::String
    neg_link_color::String
end

@option mutable struct SchematicOptions
    # relative path (within ProjectConfig.DATA_ROOT) to jld2 file containing data
    data_file::String

    # relative path (within ProjectConfig.FIG_ROOT) to save figure file
    save_file::String

    # which dataset in data_file to plot
    index::Int

    # total width of the figure (in pt)
    total_width::Float64

    # total height of the figure (in pt)
    total_height::Float64

    # list of locations to put arrowheads on switched trajectory
    # defined by fractional cumulative arc length
    arrow_locs::Vector{Float64}

    # plot energy, work, and switching signal until energy drops
    # to this point for the first time
    V_thresh::Float64

    # sizes of outer columns (panels) (relative)
    outer_colsizes::Vector{Float64}

    # column gaps between outer panels (relative)
    outer_colgaps::Vector{Float64}

    # row sizes for panel A (relative)
    panel_a_rowsizes::Vector{Float64}

    # row gap for panel A (relative)
    panel_a_rowgap::Float64

    # row sizes for panel C (relative)
    panel_c_rowsizes::Vector{Float64}

    # row gap for panel C (relative)
    panel_c_rowgap::Float64

    net::NetworkOptions
end

# c1 = colorant"rgb(223,26,54)"
# c2 = colorant"rgb(51, 80, 219)"
# r1 = range(c1, colorant"white")
# r2 = range(colorant"white", c2)
# cmap = vcat(r1, r2)
# cmap = cmap[end:-1:1]

schematic_theme = Theme(
    figure_padding = (0, 2, 0, 0),
    StreamLines = (
        resolution = (100, 100),
        density = 0.5,
        color = :black,
        linewidth = 0.5,
        arrowsize = 8,
        arrowhead = '⌃',
        minlength = 0.2,
        maxlength = 4.0,
        broken_streamlines = true
    ),
    Contourf = (
        colormap = RGBAf.(Colors.color.(to_colormap(:coolwarm)), 0.75),
        #colormap = cmap,
        levels = 20
    ),
    Lines = (
        linewidth = 0.75,
    )
)
