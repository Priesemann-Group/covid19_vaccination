# @Author: Simon Bauer
# @Date:   2021-03-16

import matplotlib as mpl

mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=["#3C627D","#49B3AC","#B0D166","#FFF576","#FB726E","#73004E"])

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

custom_cmaps = dict()
custom_cmaps["cold"] = [
    (0, "white"),
    (0.25, "#DCEDC8"),
    (0.45, "#42B3D5"),
    (0.75, "#1A237E"),
    (1, "black"),
]
custom_cmaps["hot"] = [
    (0, "white"),
    (0.3, "#FEEB65"),
    (0.65, "#cp E4521B"),
    (0.85, "#4D342F"),
    (1, "black"),
]
custom_cmaps["pinks"] = [
    (0, "white"),
    (0.2, "#FFECB3"),
    (0.45, "#E85285"),
    (0.65, "#6A1B9A"),
    (1, "black"),
]

def cmap_for_mpl(colors, n_bins=512):
    return LinearSegmentedColormap.from_list("custom_cmap", colors, N=n_bins)

cmap = cmap_for_mpl(custom_cmaps["cold"])

mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=[cmap(1/7),cmap(2/7),cmap(3/7),cmap(4/7),cmap(5/7),cmap(6/7)])