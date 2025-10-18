using Colors
using LinearAlgebra
using Makie
using Makie: MakieCore
using PyCall
using Unzip
using Infiltrator

import PyPlot

mpl_patches = pyimport("matplotlib.patches")
FancyArrowPatch = mpl_patches.FancyArrowPatch

@recipe(StreamLines, f, limits) do scene
    l_theme = default_theme(scene, Lines)
    a_theme = default_theme(scene, Arrows)
    merge(
        Attributes(
            resolution = (200, 200),
            density = 1.0,
            minlength = 0.1,
            maxlength = 4.0,
            broken_streamlines = true,

            # lines
            color = l_theme.color,
            linewidth = l_theme.linewidth,
            linestyle = l_theme.linestyle,

            # arrows
            arrowsize = a_theme.arrowsize,
            arrowhead = a_theme.arrowhead
        ),
        #default_theme(scene, Lines) # so that we can theme the lines as needed.
        #default_theme(scene, Arrows),
    )
end

function Makie.convert_arguments(::Type{<: StreamLines}, f::Function, xrange, yrange)
    xmin, xmax = extrema(xrange)
    ymin, ymax = extrema(yrange)
    return (f, Rect(xmin, ymin, xmax - xmin, ymax - ymin))
end

function Makie.convert_arguments(::Type{<: StreamLines}, f::Function, limits::Rect)
    return (f, limits)
end

function get_streamlines(f, limits::Rect{2, T}, 
                         resolution::Union{NTuple{2, Int64}, Int64},
                         density::Union{NTuple{2, Float64}, Float64},
                         minlength::Float64, maxlength::Float64,
                         broken_streamlines::Bool) where {T}
    params = zip(extrema(limits)..., Iterators.cycle(resolution))
    x_grid, y_grid = map(params) do (min_, max_, resolution_)
        LinRange(min_, max_, resolution_ + 1)
    end

    # get vector field sample for pyplot
    # XY is a Matrix (grid) of vectors [x, y], similar for UV
    XY = collect(map(collect, Iterators.product(x_grid, y_grid)))'
    UV = f.(XY)

    # dummy pyplot figure; don't need it; just want to get the streamlines
    # themselves!
    fig_dummy = PyPlot.figure()
    ax_dummy = fig_dummy.add_subplot(1, 1, 1)

    stream = ax_dummy.streamplot(first.(XY), last.(XY), first.(UV), last.(UV),
                                 density = density, minlength = minlength,
                                 maxlength = maxlength, 
                                 broken_streamlines = broken_streamlines)
    PyPlot.close(fig_dummy)

    paths = stream.lines._paths
    
    line_pts = Vector{Point2f}()

    # include the two vertices comprising the first line segment
    append!(line_pts, map(Point2f, eachrow(popfirst!(paths).vertices)))

    # @infiltrate
    # for each additional line segment...
    for path in paths
        # if end of previous line segment is not the same point as 
        # beginning of this one, insert a break (with NaNs) and 
        # include the start of the new streamline
        # if !isapprox(line_pts[end], path.vertices[1, :])
        #     push!(line_pts, Point2f(NaN, NaN))
        #     push!(line_pts, Point2f(path.vertices[1, :]))
        # end
        # # include the end pt of this line segment
        # push!(line_pts, Point2f(path.vertices[2, :]))

        # if !isapprox(line_pts[end], path.vertices[1, :])
        #     push!(line_pts, Point2f(path.vertices[1, :]))
        # end
        # include the end pt of this line segment
        append!(line_pts, Point2f.(eachrow(path.vertices)))
        push!(line_pts, Point2f(NaN, NaN))
    end

    # get FancyArrow patches
    # need to do this from the axis object, NOT stream.arrows because the latter 
    # has erased the FancyArrow information by converting them to generic Patch
    # objects
    is_fancy_arrow(p) = pybuiltin(:isinstance)(p, FancyArrowPatch)
    mpl_arrows = filter(is_fancy_arrow, ax_dummy.get_children())

    tails, tips = unzip(tuple(Point2f.(a._posA_posB)...) for a in mpl_arrows)

    # arrow_pts and arrow_dirs are each a Vector{Point2f}
    function normalized_dir(tip, tail)
        dir = f(tip)
        return 1e-6 / norm(dir) * dir
    end

    tmp = [(tip, normalized_dir(tip, tail)) for (tip, tail) in 
           zip(tips, tails) if tip != tail]
    arrow_pts, arrow_dirs = unzip(tmp)

    return (line_pts, arrow_pts, arrow_dirs)
end

function Makie.plot!(p::StreamLines)
    data = @lift begin
        get_streamlines($(p.f), $(p.limits), $(p.resolution),
                        $(p.density), $(p.minlength), $(p.maxlength),
                        $(p.broken_streamlines))
    end
    # data = lift(p.f, p.limits, p.grid_size, 
    #             p.density, p.minlength, p.maxlength) do f, limits, grid_size, density,
    #                                                     minlength, maxlength
    #     get_streamlines(f, limits, grid_size, density, minlength, maxlength)
    # end

    # function get_color(data_idx)
    #     lift(p.color_func, p.color, data) do color_func, color, data
    #         if isnothing(color_func)
    #             color
    #         else
    #             colors = color_func.(data[data_idx])
    #             colors[isnan.(data[data_idx])] .= 0.0
    #             colors
    #         end
    #     end
    # end
    
    lines_ = lines!(p, @lift($data[1]), 
                    color = p.color,
                    linestyle = p.linestyle,
                    linewidth = p.linewidth)
    

    arrows_ = arrows!(p, @lift($data[2]), @lift($data[3]),
                      color = p.color,
                      arrowsize = p.arrowsize,
                      arrowhead = p.arrowhead)
    return (lines_, arrows_)
end
