# Figure A1: Vertical gradient proxi intercomparison (MLD-sigma 26.5)

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import math
import scipy.io as sio
import matplotlib.gridspec as gridspec
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.colors as mcolors
from matplotlib.legend_handler import HandlerLine2D, HandlerTuple
from matplotlib.patches import Rectangle

# GRID:
meanm               = sio.loadmat('...\intercomparison_metrics\mean_w_womld.mat')
lonbox5             = meanm['lon_box5m']
latbox5             = meanm['lat_box5m']

# Vertical gradient:
slp                 = sio.loadmat('...\intercomparison_metrics\slope_w_55new_26_5.mat')['m_55']

# Figure:

cmap_LONM           = cm.get_cmap('BrBG_r')
boundLONM           = np.arange(-1,1.1,0.1)

boundLONM_ticks     = [-10, -8, -6, -4,-2, 0]
boundLONM_tickss    = np.multiply(boundLONM_ticks,scale)
norm_LONM           = mpl.colors.BoundaryNorm(boundLONM, cmap_LONM.N)

label_name = [r'(a) $OLIV3$',r'(b) $OMEGA3D$',r'(c) $GLORYS12v1$',r'(d) $ECCOv4r4$']

fig, axes_r = plt.subplots(nrows=4, ncols=1,sharex=True,figsize=(14, 16),subplot_kw={"projection": ccrs.Miller(central_longitude=-60)}, layout='constrained',)
proj = ccrs.Robinson(central_longitude=-60)

label_text = ['(a)','(b)','(c)','(d)']

fig.suptitle(r'Vertical gradient of $w$', fontsize = 24, y = 1.03)
ind = [0,1,3,5]
for ii,ax0 in enumerate(axes_r.flat):
    
    ax0.set_title(label_name[ii],fontsize=22)
    ax0.set_extent([-180,180,-60,60],crs=ccrs.PlateCarree())
    
    wk = slp[:,:,ind[ii]]/scale
    wk[latbox5==0] = np.nan
    wk[latbox5==-5] = np.nan
    ax0.pcolor(lonbox5+2.5,latbox5+2.5, wk, transform=ccrs.PlateCarree(), cmap = cmap_LONM, norm = norm_LONM)

# Coastlines ------------------
    land = cfeature.NaturalEarthFeature('physical', 'land', scale='50m', edgecolor='none', facecolor=cfeature.COLORS['land'], linewidth=.25)
    ax0.add_feature(land, facecolor='k')
    
# Grid ------------------------
    gl = ax0.gridlines(draw_labels=True, 
                       xlocs=range(-180, 181, 90), 
                       ylocs=range(-60, 61, 30), 
                       color='gray', zorder=1)

    gl.left_labels = True
    gl.bottom_labels = True if ii >= 3 else False     # Bottom labels for last row
    gl.top_labels = False                             # No top labels to avoid clutter 
    
    gl.xlabel_style = {"size": 16, "color": "black"}
    gl.ylabel_style = {"size": 16, "color": "black"}
    
    for spine in ax0.spines.values():
        spine.set_linewidth(2)
    ax0.tick_params(axis="both", width=2, length=6, labelsize=12)
    
# Colorbar: 
scale       = 10**(-6)

cmap_LONM           = cm.get_cmap('BrBG_r')
boundLONM           = np.arange(-1,1.1,0.1)
boundLONM_ticks     = np.arange(-1,1.2,0.2)
norm_LONM           = mpl.colors.BoundaryNorm(boundLONM, cmap_LONM.N)

cbar = plt.colorbar(
        mpl.cm.ScalarMappable(cmap=, norm = norm_LONM),
    extend='both',
    ticks=boundLONM_ticks,
    spacing='proportional',
    #orientation='horizontal',
    shrink=0.5,
    ax=axes_r,
    location='bottom',
    pad = 0.03,
)

cbar.ax.tick_params(labelsize=18)
cbar.set_label(label=r'$(10^{-6} m$ $s^{-1}$ $/$ $kg$ $m^{-3})$',size=20)

plt.savefig('.../figures/figureA1_VF.png', bbox_inches='tight', dpi=300)
