
def get_plot_layout():
    """
    Return kat's preferred plot layout options
    """

    layout = {
        "template": "simple_white",
        "margin": {
            "l": 60,
            "r": 60,
            "t": 60,
            "b": 60
        },
        "xaxis": {
            "mirror": True,
            "showgrid": True,
            "ticks": "inside",
            "automargin": False,

        },
        "yaxis": {
            "mirror": True,
            "showgrid": True,
            "ticks": "inside",
            "automargin": False,
        },
        "font": {
            "size": 16
        },
        "legend": { "xanchor": "right", "yanchor": "top", "x": 1, "y": 1 },
        "width": 600,
        "height": 600,

    }
    return layout
