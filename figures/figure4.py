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
meanm   = sio.loadmat('...\intercomparison_metrics\\mean_w.mat')
lonbox5 = meanm['lon_box5m']
latbox5 = meanm['lat_box5m']

# MLD mask 
mldmask_oli = sio.loadmat('...\OCCITENS\oliv3_glob_5_mldmask.mat')['mld_cont_bm']
mldmask_ecc = sio.loadmat('...\OCCITENS\ecco_glob_5_mldmask.mat')['mld_cont_bm'] 

# LVB mask
lvbmask_occ = sio.loadmat('...\OCCITENS\occitens_glob_5_lvbmask.mat')['lvb_bm']

# Mean:
wm    = meanm['wm']

k = 20

# Colormap:
RBensm  = sio.loadmat('figure1/RBens_new.mat')
RBens   = RBensm['RBens']

RBens1m = sio.loadmat('figure1/RBens1_new.mat')
RBens1  = RBens1m['RBens1']

cmp     = ListedColormap(RBens)
cmp1    = ListedColormap(RBens1)

# Bounds:
    
scale       = 10**(-6)
contens     = [-200, -10, -8, -6, -4, -2, -1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1, 2, 4, 6, 8, 10, 200]
contens2    = np.multiply(contens,scale)

bounds      = [-10,-8,-6,-4,-2,-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1,2,4,6,8,10]
boundticks  = [-10,-8,-6,-4,-2,-1,0,1,2,4,6,8,10]
norm        = mpl.colors.BoundaryNorm(bounds, cmp1.N)

## Quantify overstimation:
for ii in np.arange(0,2):
    for jj in np.arange(ii,6):
        diff_wm = np.nanmedian(abs(wm[:,:,k,jj])-abs(wm[:,:,k,ii]))
        print(label_name[ii],' - ',label_name[jj])
        print(diff_wm)

# Figure:

label_name = [r'(a) $OLIV3$',r'(b) $OMEGA3D$',r'(c) $GLORYS12v1$',r'(d) $ECCOv4r4$']

fig, axes_r = plt.subplots(nrows=4, ncols=1,sharex=True,figsize=(14, 16),subplot_kw={"projection": ccrs.Miller(central_longitude=-60)}, layout='constrained',)
proj = ccrs.Robinson(central_longitude=-60)

fig.suptitle(r'Time-mean $w$ at $\sigma$ 26 kg m$^{-3}$', fontsize = 24, y = 1.03)

k = 20

label_text = ['(a)','(b)','(c)','(d)']
ind = [0,1,3,5]
for ii,ax0 in enumerate(axes_r.flat):
    
    ax0.set_title(label_name[ii],fontsize=22)
    ax0.set_extent([-180,180,-60,60],crs=ccrs.PlateCarree())
    m_variable = wm[:,:,:,ind[ii]]
    
    wk = m_variable[:,:,k]
    wk[latbox5==0] = np.nan
    wk[latbox5==-5] = np.nan
    ax0.pcolor(lonbox5+2.5,latbox5+2.5,wk, transform=ccrs.PlateCarree(), cmap = cmp, vmin = -10*10**(-6), vmax = 10*10**(-6))
    
# MLD mask --------------------
    plt.rcParams['hatch.color'] = 'white'  # Choose your hatch color here
    plt.rcParams['hatch.linewidth'] = 2.0
    if ind[ii] == 2 or ind[ii] == 3:
        hatch_pcolor = ax0.pcolor(
            lonbox5 + 2.5, latbox5 + 2.5, mldmask_occ[:, :, k], 
            transform=ccrs.PlateCarree(), cmap="Grays", alpha=0, hatch="///"
        )
    if ind[ii] == 0 or ind[ii] == 1:
        hatch_pcolor = ax0.pcolor(
            lonbox5 + 2.5, latbox5 + 2.5, mldmask_oli[:, :, k], 
            transform=ccrs.PlateCarree(), cmap="Grays", alpha=0, hatch="///"
        )
    if ind[ii] == 5:
        hatch_pcolor = ax0.pcolor(
            lonbox5 + 2.5, latbox5 + 2.5, mldmask_ecc[:, :, k], 
            transform=ccrs.PlateCarree(), cmap="Grays", alpha=0, hatch="///"
        )
    plt.rcParams['hatch.color'] = 'black'
    plt.rcParams['hatch.linewidth'] = 2.0
    if ii == 0:
        lvbk = lvbmask_occ[:, :, k]
        lvbk[latbox5==0] = np.nan
        lvbk[latbox5==-5] = np.nan
        
        m_variable = wm[:,:,:,4]
        
        wk = m_variable[:,:,k]
        wk[latbox5==0] = np.nan
        wk[latbox5==-5] = np.nan
        
        lvbk[np.isnan(wk)==True] = np.nan
        hatch_pcolor = ax0.pcolor(
            lonbox5 + 2.5, latbox5 + 2.5, lvbk, 
            transform=ccrs.PlateCarree(), cmap="Grays", alpha=0, hatch="\\\\"
        )
# Coastlines ------------------
    land = cfeature.NaturalEarthFeature('physical', 'land', scale='50m', edgecolor='none', facecolor=cfeature.COLORS['land'], linewidth=.25)
    ax0.add_feature(land, facecolor='k')

# Grid ------------------------
    gl = ax0.gridlines(draw_labels=True, 
                       xlocs=range(-180, 181, 90), 
                       ylocs=range(-60, 61, 30), 
                       color='gray', zorder=1)

    gl.left_labels = True   # Left labels only for left column
    gl.bottom_labels = True if ii >= 3 else False     # Bottom labels for last row
    gl.top_labels = False                             # No top labels to avoid clutter 
    
    gl.xlabel_style = {"size": 16, "color": "black"}
    gl.ylabel_style = {"size": 16, "color": "black"}
    
    for spine in ax0.spines.values():
        spine.set_linewidth(2)
    ax0.tick_params(axis="both", width=2, length=6, labelsize=12)
    
# Colorbar 
    
scale       = 10**(-6)
contens     = [-200, -10, -8, -6, -4, -2, -1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1, 2, 4, 6, 8, 10, 200]
contens2    = np.multiply(contens,scale)

bounds      = [-10,-8,-6,-4,-2,-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1,2,4,6,8,10]
boundticks  = [-10,-8,-6,-4,-2,0,2,4,6,8,10]
norm_cont    = mpl.colors.BoundaryNorm(bounds, cmp1.N)

cbar = plt.colorbar(
        mpl.cm.ScalarMappable(cmap=cmp1, norm = norm_cont),
    extend='both',
    ticks=boundticks,
    spacing='proportional',
    #orientation='horizontal',
    shrink=0.5,
    ax=axes_r,
    location='bottom',
    pad = 0.03,
)

cbar.ax.tick_params(labelsize=18)
cbar.set_label(label=r'$10^{-6}$ $m$ $s^{-1}$',size=20)

plt.savefig('figure4_VF.png', bbox_inches='tight', dpi=300)

