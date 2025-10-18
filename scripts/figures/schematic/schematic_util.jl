using Colors
using FilePaths: PosixPath
using FilePathsBase: /
using FromFile
using Graphs
using Makie
using MakieTeX
using SimpleWeightedGraphs
using Mustache
using SparseArrays
using TemporalSync

this_dir = PosixPath(@__DIR__())

include(string(this_dir / "streamlines.jl"))
#@from "streamlines.jl" import streamlines!

@from "schematic_options.jl" import NetworkOptions

edge_label_pos_uni = Dict(
    (1, 2) => "auto=right, inner sep=1pt",
    (2, 1) => "auto=left, inner sep=1pt",
    (1, 3) => "auto=left, inner sep=1pt",
    (3, 1) => "auto=right, inner sep=1pt",
    (2, 3) => "below",
    (3, 2) => "below"
)

edge_label_pos_bi = Dict(
    (1, 2) => "auto=left, inner sep=1pt",
    (2, 1) => "auto=left, inner sep=1pt",
    (1, 3) => "auto=right, inner sep=1pt",
    (3, 1) => "auto=right, inner sep=1pt",
    (2, 3) => "below",
    (3, 2) => "above"
)

function one_edge_code(g, e)
    s = src(e)
    t = dst(e)
    loc = has_edge(g, t, s) ? edge_label_pos_bi[s, t] : edge_label_pos_uni[s, t]
    edge_style = has_edge(g, t, s) ? "auto shift" : "->"
    w = get_weight(g, s, t)
    w_trunc = trunc(Int, w)
    w = w == w_trunc ? w_trunc : w
    edge_color = w > 0 ? "positive" : "negative"
    "($s) edge[$edge_style, $edge_color] node[$loc]{$w} ($t)"
end

all_edges_code(g) =
    join(map(e -> one_edge_code(g, e), edges(g)), "\n")

function tikz_code(g, opts::NetworkOptions)
    template = read(string(this_dir / "template.tex"), String)
    d = Dict(string(key)=>getfield(opts, key) for key ∈ fieldnames(typeof(opts)))
    d["edge_code"] = all_edges_code(g)
    render(template, d, tags=("<<", ">>"))
end

net_schematic_fig(grid::GridPosition, 
                  g, 
                  opts::NetworkOptions) =
    LTeX(grid, TeXDocument(tikz_code(g, opts)), valign=:top, 
          tellwidth=false, tellheight=false)