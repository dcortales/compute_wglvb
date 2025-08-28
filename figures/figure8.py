# Figure 8: Latitudinal median correlation coefficient

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
import matplotlib.patches as mpatches
import matplotlib.colors as colors
from scipy.ndimage.interpolation import zoom
import xycmap
import pandas as pd

# GRID:
meanm   = sio.loadmat('...\intercomparison_metrics\mean_w.mat')
lonbox5 = meanm['lon_box5m']
latbox5 = meanm['lat_box5m']

wm    = meanm['wm']
    
# rho,rmse, and matsig
rho    = sio.loadmat('...\intercomparison_metrics\R_sign_RMSE_w.mat')['rho']
sig    = sio.loadmat('...\intercomparison_metrics\R_sign_RMSE_w.mat')['matsigcorr']
rmse   = sio.loadmat('...\intercomparison_metrics\R_sign_RMSE_w.mat')['RMSE']

# MLD mask 
mldmask_occ = sio.loadmat('...\OCCITENS\occitens_glob_5_mldmask.mat')['mld_cont_bm']
mldmask_oli = sio.loadmat('...\OLIV3\oliv3_glob_5_mldmask.mat')['mld_cont_bm']
mldmask_ecc = sio.loadmat('...\ECCO\ecco_glob_5_mldmask.mat')['mld_cont_bm'] 

lvbmask_occ = sio.loadmat('...\OCCITENS\occitens_glob_5_lvbmask.mat')['lvb_bm']

fig, axes_r = plt.subplots(nrows=1, ncols=1,sharex=True,figsize=(18, 14),subplot_kw={"projection": ccrs.Miller(central_longitude=-60)}, layout='constrained',)
proj = ccrs.Robinson(central_longitude=-60)
axes_r.pcolor(lonbox5+2.5,latbox5+2.5, ROCC*mask_R_OCC, transform=ccrs.PlateCarree(), cmap = cmp, norm = norm)
# %% Median r as a function of latitude
# Depth indices to show (as in original code)
depth_levels = [18, 20, 28]
depth_labels = [r'$\sigma$ 25.5 kg m$^{-3}$', r'$\sigma$ 26 kg m$^{-3}$', r'$\sigma$ 27 kg m$^{-3}$']

iii = [1,4,4,2,4,5,3,4,5,5]
jjj = [0,1,0,0,2,3,0,3,0,4]

title_name = [
    '$OLIV3$ vs $ECCO$',
    '$OLIV3$ vs $GLORYS$',
    '$OLIV3$ vs $OMEGA3D$', 
    '$OMEGA3D$ vs $ECCO$',
    '$OMEGA3D$ vs $GLORYS$',
    '$GLORYS$ vs $ECCO$'
]

linestyle = [ ':', '-',  '-','--', '--','--', ]

title_name = ['$GLORYS$ $vs$ $ECCO$',       '$OLIV3$ $vs$ $ECCO$', '$OLIV3$ $vs$ $GLORYS$',  '$OLIV3$ $vs$ $OMEGA3D$',
                     '$OMEGA3D$ $vs$ $ECCO$',
                  '$OMEGA3D$ $vs$ $GLORYS$',
              ]

iii = [5,5,3,1,5,3]
jjj = [3,0,0,0,1,1]

num_colors = 6  # number of dataset pairs
cmap = plt.cm.turbo  # choose colormap

# Sample colors evenly from the colormap range [0,1]
colors_list = [cmap(i / (num_colors - 1)) for i in range(num_colors)]

legend_handles = []
legend_labels = []

# Create the plot
fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(7, 9), sharey=True)

for idx_ax, (ax, k, label) in enumerate(zip(axes, depth_levels, depth_labels)):
    for idx, (i, j) in enumerate(zip(iii, jjj)):
        # Extract rho and compute median over longitude
        ROCC = rho[:,:,k,2,4]
        mask_R_OCC = np.full(ROCC.shape,np.nan)
        mask_R_OCC[ROCC>0.5] = 1
        rho_k = rho[:, :, k, j, i]
        rho_k[latbox5 == 0] = np.nan
        rho_k[latbox5 == -5] = np.nan
        lat_vals = latbox5[0,:]
        median_rho = np.nanmedian(rho_k*mask_R_OCC, axis=0)

        line, = ax.plot(lat_vals+2.5, median_rho, linestyle = linestyle[idx], label=title_name[idx], color=colors_list[idx], linewidth=2)
        
        if idx_ax == 0:
            legend_handles.append(line)
            legend_labels.append(title_name[idx])

    ax.grid(True)
    
    if idx_ax < 2:  # First and second panels
        ax.set_xticklabels([])
        ax.set_ylabel("")
    else:
        ax.set_xlabel('Latitude', fontsize=14)
    ax.set_xlim([-80,80])
    ax.set_ylim([-0.5,1])
    
    ax.hlines(0,-80,80,color='dimgray',linewidth = 1)
    
    ax.annotate(label,
            xy=(0.03, 0.15),  # Position in axes fraction (x=3%, y=92%)
            xycoords='axes fraction',
            fontsize=13,
            bbox=dict(boxstyle="round,pad=0.4", fc="white",alpha = 0.9, ec="lightgray", lw=1))
    
    for spine in ax.spines.values():
        spine.set_linewidth(2)
    ax.tick_params(axis='both', width=2, length=6, labelsize=12)

axes[1].set_ylabel('Median correlation coefficient', fontsize=14)

fig.legend(handles=legend_handles,
           labels=legend_labels,
           loc='lower center',
           ncol=2,
           fontsize=12,
           title="Dataset Pairs",
           title_fontsize=14,
           frameon=True,  
           fancybox=True, 
           edgecolor='gray',
           bbox_to_anchor=(0.55, -0.15))

plt.suptitle('Median correlation coefficient across various depths', fontsize=18, y=0.98)
plt.tight_layout()

plt.savefig('.../figures/figure8_VF.png', bbox_inches='tight', dpi=300)
