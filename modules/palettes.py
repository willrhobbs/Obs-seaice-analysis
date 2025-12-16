import matplotlib.colors as _mcolors
import cmocean.cm as _cm


class _palette:
    """
    Create a palette for a variable with anomaly and absolute variants.

    """

    def __init__(self, anomalies=None, absolute=None):
        if anomalies is not None:
            self.anomalies = anomalies

        if absolute is not None:
            self.absolute = absolute


def _to_colourmap(hex):
    rgb = [_mcolors.to_rgba(c) for c in hex]
    map = _mcolors.LinearSegmentedColormap.from_list("custom_cmap", rgb)
    return map


# Define hex colors
# ColourBrewer RdYlBu but without the Yl
_temperature_colours = ["#d73027", "#f46d43", "#ffffff", "#74add1", "#4575b4"]

temperature = _palette(_to_colourmap(_temperature_colours), absolute=_cm.thermal)


# ColourBrewer BrBG
_sailnity_anomalies = [
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


sailnity = _palette(anomalies=_to_colourmap(_sailnity_anomalies), absolute=_cm.haline)

# _sea_ice_absolute = ["#040613", "#3c4185", "#4177b6", "#7ec1d1", "#f4ffff"]
# Use transparency
_sea_ice_absolute = ["#FFFFFF00", "#FFFFFFFF"]

_sea_ice_anomalies = [
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

sea_ice = _palette(
    anomalies=_to_colourmap(_sea_ice_anomalies),
    absolute=_to_colourmap(_sea_ice_absolute),
)
