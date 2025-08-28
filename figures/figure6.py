# Figure 6: Variance intercomparison
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
meanm   = sio.loadmat('...\intercomparison_metrics\mean_w.mat')
lonbox5 = meanm['lon_box5m']
latbox5 = meanm['lat_box5m']
    
# Variance:
wv      = np.squeeze(sio.loadmat('...\intercomparison_metrics\variance_w.mat')['wv'])

# MLD mask 
mldmask_occ = sio.loadmat('...\OCCITENS\occitens_glob_5_mldmask.mat')['mld_cont_bm']
mldmask_oli = sio.loadmat('...\OLIV3\oliv3_glob_5_mldmask.mat')['mld_cont_bm']
mldmask_ecc = sio.loadmat('...\ECCO\ecco_glob_5_mldmask.mat')['mld_cont_bm'] 

lvbmask_occ = sio.loadmat('...\OCCITENS\occitens_glob_5_lvbmask.mat')['lvb_bm']

# Correlation coefficient
rho    = sio.loadmat('...\intercomparison_metrics\R_sign_RMSE_w.mat')['rho']

k = 20

label_name = [r'$OLIV3$',r'$OMEGA3D$',r'$NEMO$',r'$GLORYS12v1$',r'$NEMO$ ($w_g$)',r'$ECCOv4r4$']

# Figure:

cmap3           = cm.get_cmap('plasma')
bounds_var      = [-14,-13.8,-13.6,-13.4,-13.2,-13,-12.8,-12.6,-12.4,-12.2,-12,-11.8,-11.6,-11.4,-11.2,-11]
norm2 = mpl.colors.BoundaryNorm(bounds_var, cmap3.N)
boundticks_var  = [-14,-13,-12,-11]

label_name = [r'(a) $OLIV3$',r'(b) $OMEGA3D$',r'(c) $GLORYS12v1$',r'(d) $ECCOv4r4$']

fig, axes_r = plt.subplots(nrows=4, ncols=1,sharex=True,figsize=(14, 16),subplot_kw={"projection": ccrs.Miller(central_longitude=-60)}, layout='constrained',)
proj = ccrs.Robinson(central_longitude=-60)

fig.suptitle(r'$w$ variance at $\sigma$ 26 kg m$^{-3}$', fontsize = 24, y = 1.03)

k = 20

label_text = ['(a)','(b)','(c)','(d)']
ind = [0,1,3,5]
for ii,ax0 in enumerate(axes_r.flat):
    
    ax0.set_title(label_name[ii],fontsize=22)
    ax0.set_extent([-180,180,-60,60],crs=ccrs.PlateCarree())
    m_variable = wv[:,:,:,ind[ii]]
    
    wk = m_variable[:,:,k]
    wk[latbox5==0] = np.nan
    wk[latbox5==-5] = np.nan
    ax0.pcolor(lonbox5+2.5,latbox5+2.5,wk, transform=ccrs.PlateCarree(), cmap = 'plasma', norm = norm2)

    plt.rcParams['hatch.color'] = 'black'
    plt.rcParams['hatch.linewidth'] = 2.0
    if ii == 0 or ii == 4:
        lvbk = lvbmask_occ[:, :, k]
        lvbk[latbox5==0] = np.nan
        lvbk[latbox5==-5] = np.nan
        
        m_variable = wv[:,:,:,4]
        
        wk = m_variable[:,:,k]
        wk[latbox5==0] = np.nan
        wk[latbox5==-5] = np.nan
        
        r_mask = np.full(wk.shape, np.nan)
        r_mask[rho[:,:,k,2,4] < 0.5] = 1 # Mask correlation coefficient nemo vs nemo wg
        r_mask[latbox5==0] = np.nan
        r_mask[latbox5==-5] = np.nan
        hatch_pcolor = ax0.pcolor(
            lonbox5 + 2.5, latbox5 + 2.5, r_mask, 
            transform=ccrs.PlateCarree(), cmap="Grays", alpha=0, hatch="\\\\"
        )

# Coastlines ------------------
    land = cfeature.NaturalEarthFeature('physical', 'land', scale='50m', edgecolor='none', facecolor=cfeature.COLORS['land'], linewidth=.25)
    ax0.add_feature(land, facecolor='k')
    
# Grid ------------------------
    gl = ax0.gridlines(draw_labels=True, xlocs=range(-180, 181, 90), ylocs=range(-60, 61, 30), color='gray', zorder=1)

    gl.left_labels = True 
    gl.bottom_labels = True if ii >= 4 else False   
    gl.top_labels = False                             
    
    gl.xlabel_style = {"size": 16, "color": "black"}
    gl.ylabel_style = {"size": 16, "color": "black"}
    
    for spine in ax0.spines.values():
        spine.set_linewidth(2)
    ax0.tick_params(axis="both", width=2, length=6, labelsize=12)

# Colorbar:

cbar = plt.colorbar(
        mpl.cm.ScalarMappable(cmap='plasma',norm = norm2),
    extend='both',
    ticks=boundticks_var,
    spacing='proportional',
    shrink=0.5,
    ax=axes_r,
    location='bottom',
    pad = 0.03,
)

cbar.ax.tick_params(labelsize=18)
cbar.set_label(label=r'$variance$ $[log_{10}(m^2$ $s^{-2)}]$',size=20)


plt.savefig('figure6_VF.png', bbox_inches='tight', dpi=300)

