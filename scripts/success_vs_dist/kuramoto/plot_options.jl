using Configurations

@option struct PanelAOptions
    tick_pad::Float64
    xticks::Vector{Float64}
    panel_label_x::Float64
    marker::String
    marker_size::Float64
    x_axis_label_pad::Float64
end

@option struct LegendOptions
    marker::String
    marker_size::Float64
    loc_x::Float64
    loc_y::Float64
    title_x::Float64
    title_y::Float64
    title_size::Float64
end

@option struct PanelBOptions
    elev::Float64
    azim::Float64
    dist::Float64
    mode::String
    alpha::Float64
    tick_pad::Float64
    bins::Int64
    bandwidth::Float64
    x_axis_label_pad::Float64
    z_axis_label_pad::Float64
    legend_label_spacing::Float64
    handle_text_pad::Float64

    panel_label_x::Float64
end

@option struct PlotOptions
    V₀_bins::Int64
    cmap::String
    m::Int64
    n::Int64
    p_range::Vector{Float64}
    fig_file_name::String
    width_ratios::Vector{Float64}

    fig_width::Float64
    fig_height::Float64
    tick_label_size::Float64
    major_tick_size::Float64
    minor_tick_size::Float64
    font_size::Float64
    axis_label_pad::Float64
    axis_label_size::Float64

    panel_label_size::Float64
    panel_label_y::Float64

    panelA::PanelAOptions
    panelB::PanelBOptions
    legend::LegendOptions
end