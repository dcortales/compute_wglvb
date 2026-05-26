# FIGURE 5: Time-mean intercomparison

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.io as sio
from matplotlib.colors import ListedColormap
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# GRID:
meanm   = sio.loadmat('...\intercomparison_metrics\mean_w.mat')
lonbox5 = meanm['lon_box5m']
latbox5 = meanm['lat_box5m']

# MLD mask 
mldmask_occ = sio.loadmat('...\Data\\OGCM\\occitens_glob_5_mldmask.mat')['mld_cont_bm']
mldmask_oli = sio.loadmat('...\\Data\\OGCM\\oliv3_glob_5_mldmask.mat')['mld_cont_bm']
mldmask_ecc = sio.loadmat('...\\Data\\OGCM\\ecco_glob_5_mldmask.mat')['mld_cont_bm'] 

lvbmask_occ = sio.loadmat('...\\Data\\OGCM\\occitens_glob_5_lvbmask.mat')['lvb_bm']

# Mean:
wm    = meanm['wm']

k = 20

# Colormap:
RBens1m = sio.loadmat('figure1/RBens1_new.mat')
RBens1  = RBens1m['RBens1']
cmp1    = ListedColormap(RBens1)

# Bounds:
scale = 0.1
bounds      = [-10,-8,-6,-4,-2,-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1,2,4,6,8,10]
bounds2     = np.multiply(bounds,scale)
norm        = mpl.colors.BoundaryNorm(bounds2, cmp1.N)

# %% Figure
label_name = [r'(a) $OLIV3$',r'(b) $OMEGA3D$',r'(c) $GLORYS12v1$',r'(d) $ECCOv4r4$',
              r'(e) |$OMEGA3D$| - $|OLIV3|$',r'(f) $|GLORYS12v1|$ - $|OLIV3|$',r'(g) $|ECCOv4r4|$ - $|OLIV3|$',
              r'(h) $|GLORYS12v1|$ - $|OMEGA3D|$',r'(i) $|ECCOv4r4|$ - $|OMEGA3D|$', r'(j) $|ECCOv4r4|$ - $|GLORYS12v1|$']

fig, axes_r = plt.subplots(nrows=5, ncols=2,sharex=True,figsize=(16, 18),subplot_kw={"projection": ccrs.Miller(central_longitude=-60)}, layout='constrained',)
proj = ccrs.Robinson(central_longitude=-60)

fig.suptitle(r'Time-mean $w$ and differences among products at $\sigma$ 26 kg m$^{-3}$', fontsize = 24, y = 1.03)

k = 20

ind = [0,1,3,5,1,3,5,3,5,5]
for ii,ax0 in enumerate(axes_r.flat):
    
    ax0.set_title(label_name[ii],fontsize=22)
    ax0.set_extent([-180,180,-60,60],crs=ccrs.PlateCarree())
    if ii < 4:
        m_variable = wm[:,:,:,ind[ii]]
    elif ii > 3 and ii < 7:
        m_variable = abs(wm[:,:,:,ind[ii]]) - abs(wm[:,:,:,0])
    elif ii > 6 and ii < 9:
        m_variable = abs(wm[:,:,:,ind[ii]]) - abs(wm[:,:,:,1])
    elif ii > 8:
        m_variable = abs(wm[:,:,:,ind[ii]]) - abs(wm[:,:,:,3])
    
    wk = m_variable[:,:,k]
    wk[latbox5==0] = np.nan
    wk[latbox5==-5] = np.nan
    ax0.pcolor(lonbox5+2.5,latbox5+2.5,wk*60*60*24, transform=ccrs.PlateCarree(), cmap = cmp1, norm = norm)

# MLD mask --------------------
    if ii < 4:
        plt.rcParams['hatch.color'] = 'white' 
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
    gl.left_labels = True 
    gl.bottom_labels = True if ii >= 7 else False  
    gl.top_labels = False                     
    
    gl.xlabel_style = {"size": 16, "color": "black"}
    gl.ylabel_style = {"size": 16, "color": "black"}
    
    for spine in ax0.spines.values():
        spine.set_linewidth(2)
    ax0.tick_params(axis="both", width=2, length=6, labelsize=12)
    
# Colorbar 
scale = 0.1
bounds       = [-10,-8,-6,-4,-2,-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1,2,4,6,8,10]
bounds2      = np.multiply(bounds,scale)
bounds       = [-10,-8,-6,-4,-2,0,2,4,6,8,10]
boundticks2  = np.multiply(bounds,scale)
norm_cont    = mpl.colors.BoundaryNorm(bounds2, cmp1.N)

cbar = plt.colorbar(
        mpl.cm.ScalarMappable(cmap=cmp1, norm = norm_cont),
    extend='both',
    ticks=boundticks2,
    spacing='proportional',
    shrink=0.5,
    ax=axes_r,
    location='bottom',
    pad = 0.03,
)

cbar.ax.tick_params(labelsize=18)
cbar.set_label(label=r'$m$ $day^{-1}$',size=20)

plt.savefig('...\figures\figure5.png', bbox_inches='tight', dpi=300)
