#!/usr/bin/env python
"""GLobally used variables."""
from matplotlib import colors
import numpy as np

r=np.linspace(0,0.14,10)
g=np.linspace(0,1,10)
b=np.linspace(0,0,10)
cmap_gfp_list=[l for l in zip(r,g,b)]
cmap_gfp = colors.LinearSegmentedColormap.from_list("gfp",cmap_gfp_list)