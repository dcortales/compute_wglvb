import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.io as sio
from matplotlib import cm
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.colors as colors

# GRID:
meanm   = sio.loadmat('C:\\Users\\yago_\\Documents\\LOCEAN\\OLIV3 paper\\V4 codes\\intercomparison_metrics\\mean_w.mat')
lonbox5 = meanm['lon_box5m']
latbox5 = meanm['lat_box5m']
    
# Variance:
wv      = np.squeeze(sio.loadmat('C:\\Users\\yago_\\Documents\\LOCEAN\\OLIV3 paper\\V4 codes\\intercomparison_metrics\\variance_w.mat')['wv'])

wv      = np.squeeze(sio.loadmat('C:/Users/yago_/Documents/LOCEAN/intercomparison_metrics/variance_w_nolog.mat')['wv'])
# MLD mask 
mldmask_occ = sio.loadmat('C:\\Users\\yago_\\Documents\\LOCEAN\\Data\\OCCITENS\\occitens_glob_5_mldmask.mat')['mld_cont_bm']
mldmask_oli = sio.loadmat('C:\\Users\\yago_\\Documents\\LOCEAN\\Data\\OCCITENS\\oliv3_glob_5_mldmask.mat')['mld_cont_bm']
mldmask_ecc = sio.loadmat('C:\\Users\\yago_\\Documents\\LOCEAN\\Data\\OCCITENS\\ecco_glob_5_mldmask.mat')['mld_cont_bm'] 

lvbmask_occ = sio.loadmat('C:\\Users\\yago_\\Documents\\LOCEAN\\Data\\OCCITENS\\occitens_glob_5_lvbmask.mat')['lvb_bm']

rho    = sio.loadmat('C:\\Users\\yago_\\Documents\\LOCEAN\\OLIV3 paper\\V4 codes\\intercomparison_metrics\\R_sign_RMSE_w.mat')['rho']

k = 20
# %% Functions

def truncate_cmap(cmap, minval=0.05, maxval=0.95, n=256):
    """Return a truncated copy of a colormap."""
    new_colors = cmap(
        np.linspace(minval, maxval, n)
    )
    return mpl.colors.LinearSegmentedColormap.from_list(
        f"trunc({cmap.name},{minval},{maxval})", new_colors
    )

# %% Figure 7

# Colormap
bounds      = np.concatenate([np.linspace(0, 0.01, 21),np.linspace(0.01, 0.05, 5)[1:]])
bound_names = np.concatenate([np.linspace(0, 0.01, 2), np.linspace(0.01, 0.05, 5)[1:]])
orig_cmap   = mpl.colormaps.get_cmap('gist_earth')
cmap        = truncate_cmap(orig_cmap, 0.01, 0.95)
norm        = colors.BoundaryNorm(boundaries=bounds, ncolors=cmap.N)

# Plot Figure
label_name  = [r'(a) $OLIV3$',r'(b) $OMEGA3D$',r'(c) $GLORYS12v1$',r'(d) $ECCOv4r4$']
label_text  = ['(a)','(b)','(c)','(d)']
k           = 20
ind         = [0,1,3,5]

fig, axes_r = plt.subplots(nrows=4, ncols=1,sharex=True,figsize=(14, 16),subplot_kw={"projection": ccrs.Miller(central_longitude=-60)}, layout='constrained',)
proj = ccrs.Robinson(central_longitude=-60)

fig.suptitle(r'$w$ variance at $\sigma$ 26 kg m$^{-3}$', fontsize = 24, y = 1.03)

for ii,ax0 in enumerate(axes_r.flat):
    
    ax0.set_title(label_name[ii],fontsize=22)
    ax0.set_extent([-180,180,-60,60],crs=ccrs.PlateCarree())
    m_variable = wv[:,:,:,ind[ii]]
    
    wk = m_variable[:,:,k]
    wk[latbox5==0] = np.nan
    wk[latbox5==-5] = np.nan
    ax0.pcolor(lonbox5+2.5,latbox5+2.5,wk, transform=ccrs.PlateCarree(), cmap = cmap, norm = norm)

    plt.rcParams['hatch.color'] = 'lightgray'
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
        r_mask[rho[:,:,k,2,4] < 0.5] = 1
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
    gl = ax0.gridlines(draw_labels=True, 
                       xlocs=range(-180, 181, 90), 
                       ylocs=range(-60, 61, 30), 
                       color='gray', zorder=1)

    gl.left_labels      = True  
    gl.bottom_labels    = True if ii >= 4 else False   
    gl.top_labels       = False                        
    gl.xlabel_style     = {"size": 16, "color": "black"}
    gl.ylabel_style     = {"size": 16, "color": "black"}
    
    for spine in ax0.spines.values():
        spine.set_linewidth(2)
    ax0.tick_params(axis="both", width=2, length=6, labelsize=12)

# Colorbar:
cbar = plt.colorbar(
        mpl.cm.ScalarMappable(cmap=cmap, norm=norm),
    extend='max',
    ticks=bound_names,
    spacing='proportional',
    shrink=0.5,
    ax=axes_r,
    location='bottom',
    pad = 0.03,
)

cbar.ax.tick_params(labelsize=18)
cbar.set_label(label=r'$m^2$ $day^{-2}$',size=20)

#plt.savefig('C:\\Users\\yago_\\Documents\\LOCEAN\\OLIV3 paper\\V4 codes\\review\\6_variance_w_55_VF.png', bbox_inches='tight', dpi=300)


# %% Figure 8

# Colormap
cmap3           = cm.get_cmap('PuOr_r')
bounds_var      = np.arange(-0.01,0.011,0.001)
norm2 = mpl.colors.BoundaryNorm(bounds_var, cmap3.N)
boundticks_var  = np.arange(-0.01,0.014,0.004)

# Figure
label_name  = [r'(a) $OMEGA3D$ - $OLIV3$',r'(b) $GLORYS12v1$ - $OLIV3$',r'(c) $ECCOv4r4$ - $OLIV3$',
              r'(d) $GLORYS12v1$ - $OMEGA3D$',r'(e) $ECCOv4r4$ - $OMEGA3D$', r'(f) $GLORYS12v1$ - $ECCOv4r4$']
k           = 20
ind         = [1,3,5,3,5,5]

fig, axes_r = plt.subplots(nrows=3, ncols=2,sharex=True,figsize=(14, 12),subplot_kw={"projection": ccrs.Miller(central_longitude=-60)}, layout='constrained',)
proj = ccrs.Robinson(central_longitude=-60)

fig.suptitle(r'$w$ variance differences at $\sigma$ 26 kg m$^{-3}$', fontsize = 24, y = 1.03)

for ii,ax0 in enumerate(axes_r.flat):
    
    ax0.set_title(label_name[ii],fontsize=22)
    ax0.set_extent([-180,180,-60,60],crs=ccrs.PlateCarree())
    m_variable = wv[:,:,:,ind[ii]]
    
    if  ii < 3:
        m_variable = wv[:,:,:,ind[ii]] - wv[:,:,:,0]
    elif ii > 2 and ii < 5:
        m_variable = wv[:,:,:,ind[ii]] - wv[:,:,:,1]
    elif ii > 4:
        m_variable = wv[:,:,:,ind[ii]] - wv[:,:,:,3]
    
    wk = m_variable[:,:,k]
    wk[latbox5==0] = np.nan
    wk[latbox5==-5] = np.nan
    ax0.pcolor(lonbox5+2.5,latbox5+2.5,wk, transform=ccrs.PlateCarree(), cmap = cmap3, norm = norm2)

# Coastlines ------------------
    land = cfeature.NaturalEarthFeature('physical', 'land', scale='50m', edgecolor='none', facecolor=cfeature.COLORS['land'], linewidth=.25)
    ax0.add_feature(land, facecolor='k')

# Grid ------------------------
    gl = ax0.gridlines(draw_labels=True, 
                       xlocs=range(-180, 181, 90), 
                       ylocs=range(-60, 61, 30), 
                       color='gray', zorder=1)
    gl.left_labels      = True  
    gl.bottom_labels    = True if ii >= 4 else False   
    gl.top_labels       = False                            
    gl.xlabel_style     = {"size": 16, "color": "black"}
    gl.ylabel_style     = {"size": 16, "color": "black"}
    
    for spine in ax0.spines.values():
        spine.set_linewidth(2)
    ax0.tick_params(axis="both", width=2, length=6, labelsize=12)

# Colorbar:

cbar = plt.colorbar(
        mpl.cm.ScalarMappable(cmap=cmap3,norm = norm2),
    extend='both',
    ticks=boundticks_var,
    spacing='proportional',
    shrink=0.5,
    ax=axes_r,
    location='bottom',
    pad = 0.03,
)

cbar.ax.tick_params(labelsize=18)
cbar.set_label(label=r'$m^2$ $day^{-2}$',size=20)

#plt.savefig('C:\\Users\\yago_\\Documents\\LOCEAN\\OLIV3 paper\\V4 codes\\review\\6_variance_w_55_VF_diff.png', bbox_inches='tight', dpi=300)
