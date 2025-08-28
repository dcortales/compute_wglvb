# Figure 7: Correlation cofficient intercomparison
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

# Figure:

fig, axes_r = plt.subplots(nrows=3, ncols=2,sharex=True,figsize=(18, 14),subplot_kw={"projection": ccrs.Miller(central_longitude=-60)}, layout='constrained',)
proj = ccrs.Robinson(central_longitude=-60)

fig.suptitle('$w$ correlation coefficient at $\sigma$ 26 kg m$^{-3}$', fontsize = 28, y = 1.03)

k = 20

label_text = ['(a)','(b)','(c)','(d)','(e)','(f)']
title_name = ['(a) $GLORYS$ $vs$ $ECCO$',       '(b) $OLIV3$ $vs$ $OMEGA3D$',
              '(c) $OLIV3$ $vs$ $ECCO$',        '(d) $OMEGA3D$ $vs$ $ECCO$',
              '(e) $OLIV3$ $vs$ $GLORYS$',      '(f) $OMEGA3D$ $vs$ $GLORYS$',
              ]

iii = [5,1,5,5,3,3]
jjj = [3,0,0,1,0,1]

cont_rho    = [-1, 0, 0.5,0.6,0.7,0.8,0.9,1]
boundstick  = np.arange(-1,1.2,0.2)
bounds      = [-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]
norm        = mpl.colors.BoundaryNorm(bounds, cmp.N)

for ij, ax0 in enumerate(axes_r.flat):
    
    # Field selection: 
        
        ax0.set_title(title_name[ij],fontsize=26)
        ax0.set_xlabel(r'$Longitude(^o)$',fontsize=24)
        ax0.set_ylabel(r'$Latitude(^o)$',fontsize=24)
        title_text = title_name[ij]
        Rvar, RMSEvar = rho[:,:,k,jjj[ij],iii[ij]], rmse[:,:,k,iii[ij],jjj[ij]]*10**(7)
        
        Rvar[latbox5==0] = np.nan
        Rvar[latbox5==-5] = np.nan
        
        RMSEvar[latbox5==0] = np.nan
        RMSEvar[latbox5==-5] = np.nan
        
        ax0.set_extent([-180,180,-60,60],crs=ccrs.PlateCarree())
        
        ax0.pcolor(lonbox5+2.5,latbox5+2.5, Rvar, transform=ccrs.PlateCarree(), cmap = cmp, norm = norm)
        
      
        sigk = sig[:,:,k,jjj[ij],iii[ij]]
        sigk[latbox5==0] = np.nan
        sigk[latbox5==-5] = np.nan
        ax0.scatter(lonbox5*sigk+2.5,latbox5*sigk+2.5, s=5, c='black', marker='.', transform=ccrs.PlateCarree(),zorder=1000)

    # MLD mask --------------------
        plt.rcParams['hatch.color'] = 'black'  # Choose your hatch color here
        plt.rcParams['hatch.linewidth'] = 2.0
        r_mask = np.full(sigk.shape, np.nan)
        r_mask[rho[:,:,k,2,4] < 0.5] = 1 # Mask correlation coefficient nemo vs nemo wg
        r_mask[latbox5==0] = np.nan
        r_mask[latbox5==-5] = np.nan
        if ij == 1 or ij == 2 or ij == 4:
            hatch_pcolor = ax0.pcolor(
                lonbox5 + 2.5, latbox5 + 2.5, r_mask, 
                transform=ccrs.PlateCarree(), cmap="Grays", alpha=0, hatch="\\\\"
            )

    # Coastlines ------------------
        land = cfeature.NaturalEarthFeature('physical', 'land', scale='50m', edgecolor='none', facecolor=cfeature.COLORS['land'], linewidth=.25)
        ax0.add_feature(land, facecolor='k',zorder=10000)

    # Grid ------------------------
        gl = ax0.gridlines(draw_labels=True, 
                           xlocs=range(-180, 181, 90), 
                           ylocs=range(-60, 61, 30), 
                           color='gray', zorder=1)

        gl.right_labels = False if ij % 2 == 0 else True  # Right labels only for right column
        gl.left_labels = True if ij % 2 == 0 else False   # Left labels only for left column
        gl.bottom_labels = True if ij >= 4 else False     # Bottom labels for last row
        gl.top_labels = False                             # No top labels to avoid clutter 
        
        gl.xlabel_style = {"size": 18, "color": "black"}
        gl.ylabel_style = {"size": 18, "color": "black"}
        
        for spine in ax0.spines.values():
            spine.set_linewidth(2)
        ax0.tick_params(axis="both", width=2, length=6, labelsize=12)

        # Colorbar:
cbar = plt.colorbar(
        mpl.cm.ScalarMappable(norm = norm,cmap=cmp),
    extend='both',
    ticks=boundstick,
    spacing='proportional',
    shrink=0.5,
    ax=axes_r,
    location='bottom',
)

cbar.ax.tick_params(labelsize=20)
cbar.set_label('Correlation coefficient', fontsize=22)
plt.savefig('.../figures/figure7_VF.png', bbox_inches='tight', dpi=300)

