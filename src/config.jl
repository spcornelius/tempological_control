module ProjectConfig
    using Colors
    using FilePaths: PosixPath
    using FilePathsBase: /
    using Makie

    const ALMOST_BLACK = colorant"#262626"
    const PROJECT_ROOT = PosixPath(pkgdir(@__MODULE__))
    const FIG_ROOT = PROJECT_ROOT / "figures"
    const DATA_ROOT = PROJECT_ROOT / "data"

    base_theme = Theme(
        font = "Arial", 
        fontsize = 10,
        fontcolor = ALMOST_BLACK,
        textcolor = ALMOST_BLACK,
        figure_padding = (0, 0, 0, 0),
        backgroundcolor = :white,
        Axis = (
            titlesize = 12,
            xlabelsize = 12,
            ylabelsize = 12,
            topspinecolor = ALMOST_BLACK,
            bottomspinecolor = ALMOST_BLACK,
            rightspinecolor = ALMOST_BLACK,
            leftspinecolor = ALMOST_BLACK,
            xtickcolor = ALMOST_BLACK,
            ytickcolor = ALMOST_BLACK,
            yticklabelcolor = ALMOST_BLACK,
            xticklabelcolor = ALMOST_BLACK,
            xgridvisible = false,
            ygridvisible = false,
            xticksize = 2,
            yticksize = 2,
            xtickwidth = 0.5,
            ytickwidth = 0.5,
            xgridwidth = 0.5,
            ygridwidth = 0.5,
            spinewidth = 0.5,
            xlabelpadding = 2,
            ylabelpadding = 15
        ),
        Arrows = (
            arrowhead='⌃',
        ),
        Label = (
            font = "Arial",
            fontsize = 14,
            color = ALMOST_BLACK,
            halign = :left
        ),
        Legend = (
            framevisible=false,
        )
    )
end
