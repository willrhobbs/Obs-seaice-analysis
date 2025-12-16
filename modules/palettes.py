import matplotlib.colors as mcolors
import cmocean


class palette:
    """
    Create a palette for a variable with anomaly and absolute variants.

    """

    def __init__(self, anomalies=None, absolute=None):
        if anomalies is not None:
            self.anomalies = anomalies

        if absolute is not None:
            self.absolute = absolute


def to_colourmap(hex):
    rgb = [mcolors.to_rgb(c) for c in hex]
    map = mcolors.LinearSegmentedColormap.from_list("custom_cmap", rgb)
    return map


# Define hex colors
# ColourBrewer RdYlBu but without the Yl
temperature_colours = ["#d73027", "#f46d43", "#ffffff", "#74add1", "#4575b4"]
temperature = palette(to_colourmap(temperature_colours), absolute=cmocean.cm.thermal)


temperature.absolute.__str__
# ColourBrewer BrBG
sailnity_anomalies = [
    "#8c510a",
    "#bf812d",
    "#dfc27d",
    "#f6e8c3",
    "#f5f5f5",
    "#c7eae5",
    "#80cdc1",
    "#35978f",
    "#01665e",
]


sailnity = palette(
    anomalies=to_colourmap(sailnity_anomalies), absolute=cmocean.cm.haline
)

# From the PolarWatch scale (just colourpicking the colourbar)
sea_ice_absolute = ["#040613", "#3c4185", "#4177b6", "#7ec1d1", "#f4ffff"]

sea_ice_anomalies = [
    "#b35806",
    "#e08214",
    "#fdb863",
    "#fee0b6",
    "#f7f7f7",
    "#d8daeb",
    "#b2abd2",
    "#8073ac",
    "#542788",
]
sea_ice = palette(
    anomalies=to_colourmap(sea_ice_anomalies), absolute=to_colourmap(sea_ice_absolute)
)
