using Colors
using FromFile
using Graphs
using Makie
using Makie: MakieCore
using PyCall
using SimpleWeightedGraphs
using SparseArrays
using TemporalSync
using Unzip

@from "streamlines.jl" import StreamLines, streamlines!

@recipe(VectorField, g) do scene
    s_theme = default_theme(scene, StreamLines)
    merge(
        Attributes(
            colorscale = Makie.automatic,
            resolution = s_theme.resolution,
        ),
        default_theme(scene, StreamLines), # so that we can theme the lines as needed.
        default_theme(scene, Contourf)
    )
end

meshgrid(x, y) = collect(map(collect, Iterators.product(x, y)))'

function Makie.plot!(p::VectorField)
    sys = lift(p, p.g) do g
        KuramotoSys(sparse(weights(g)))
    end

    V = KuramotoEnergy()

    ϕ_min = -pi + 0.001
    ϕ_max = pi - 0.001
    ϕ_range = ϕ_min..ϕ_max
    ϕ_grid = lift(p, p.resolution) do resolution
        resolution = to_ndim(Vec{2, Int}, resolution, last(resolution))
        meshgrid(LinRange(ϕ_min, ϕ_max, resolution[1]), 
                 LinRange(ϕ_min, ϕ_max, resolution[2]))
    end

    # project a point in 3D phase space to 2D by measuring
    # the first two angles relative to the third, i.e.
    # θ = [θ₁, θ₂, θ₃] -> ϕ = [θ₁ - θ₃, θ₂ - θ₃]
    project2d(θ) = (θ .- θ[3])[1:2]

    # reverse of the above; turn a 2D point in the projected
    # space to a 3D vector (by appending θ₃ = 0)
    to_3d(ϕ) = [ϕ..., 0.0]

    # sys is 3D system, so the dynamics (rhs), energy, and work
    # expect a 3D vector θ. Wrap them so we can calcuale on a
    # 2D point ϕ

    E(ϕ) = energy(V, to_3d(ϕ))
    f = @lift begin
        ϕ -> Point2f(project2d(rhs(to_3d(ϕ), $sys, nothing))...)
    end

    W = @lift begin
        ϕ -> work(V, to_3d(ϕ), $sys, nothing)
    end

    W_grid = @lift($W.($ϕ_grid))

    W_mag = @lift(maximum(abs.($W_grid)))
    #colorrange = 
    #W_grid, W_mag, colorrange
    # cmap = RGBAf.(Colors.color.(to_colormap(:coolwarm)), 0.33)
    x = @lift(vec(first.($ϕ_grid)))
    y = @lift(vec(last.($ϕ_grid)))
    z = @lift(vec($W_grid))

    replace_automatic!(p, :colorrange) do 
        @lift((-$W_mag, $W_mag))
    end

    # c = heatmap!(p, ϕ_range, ϕ_range, W_grid, colorrange=colorrange)
    c = contourf!(p, x, y, z, colorscale = p.colorrange)
    #c = contourf!(p, x, y, z)


    # s = streamlines!(ax, f, ϕ₁_range, ϕ₂_range, color=:black, grid_size=(200, 200), 
    #                  density=opts.density[i], linewidth=opts.linewidth, 
    #                  arrowsize=opts.arrowsize, minlength=opts.minlength, maxlength=opts.maxlength,
    #                  broken_streamlines=opts.broken_streamlines)
    s = streamlines!(p, f, ϕ_range, ϕ_range)              
    p
end