# FIGURE 1: Time-mean OLIV3 across isopycnal levels
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
import scipy.io as sio
from matplotlib.colors import ListedColormap
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import xarray as xr

# %% Load variables
# MLD
MLDm            = sio.loadmat('...\MLD_ARMOR3D_GLOB.mat')
mldcont         = MLDm['mld_contac']

# Isopycnal level depth:
filo            = '...\Data\OLIV3\oliv3_glob_025_ekf_annual_isolevm.nc'
file            = xr.open_dataset(filo)
h_isolev        = file.variables['h_isolev'].values.T
h_isolev        = h_isolev[:,:,:,0]

# OLIV3 vertical velocities and 3d grid
filo            = '...\Data\OLIV3\oliv3_glob_025_ekf_annual_isolevm_5filtr.nc'
file            = xr.open_dataset(filo)
lon, lat, iso   = file.variables['longitude'].values.T, file.variables['latitude'].values.T, file.variables['isolev'].values.T
w               = np.mean(file.variables['w_oliv3_isolevf'].values.T,3)

dimlon, dimlat, dimdep = lon.shape[0], lat.shape[1], len(iso)

circshift_lim   = int(len(lon)/2)

# Colormap definition
RBens1m     = sio.loadmat('figure1/RBens1.mat')
RBens1      = RBens1m['RBens1']
cmp1        = ListedColormap(RBens1)

# %% Figure
scale       = 0.1
bounds      = [-10,-8,-6,-4,-2,-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1,2,4,6,8,10]
bounds2     = np.multiply(bounds,scale)
norm        = mpl.colors.BoundaryNorm(bounds2, cmp1.N)

proj        = ccrs.Robinson(central_longitude=-60)

fig, axs    = plt.subplots(4, 2, figsize=(17, 19),subplot_kw={"projection": ccrs.Robinson(central_longitude=-60)}, layout='constrained',)

# oliv3:
k = [14,16,18,20,22,25,28,32]
ll = ['(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)']
for ii,ax0 in enumerate(axs.flat):
    
    ax0.set_title(ll[ii]+' $\sigma$ '+str(round(iso[k[ii]],2))+ ' kg m$^{-3}$',fontsize = 24)
    ax0.set_xlabel(r'$Longitude(^o)$')
    ax0.set_ylabel(r'$Latitude(^o)$')
    
    grd = ax0.gridlines(
        draw_labels=True, 
        xlocs=range(-180, 181, 90), 
        ylocs=range(-60, 61, 30), 
        color='gray',
        zorder=10000,
    )
    grd.top_labels = False
    grd.left_labels = False
    grd.ylabel_style = {'color': 'white'}
    grd.bottom_labels = False
    grd.top_labels = False
    
    wk = w[:,:,k[ii]]
    wk[abs(lat)<5] = np.nan
    cs = ax0.contourf(lon, lat, wk*60*60*24,bounds2,transform=ccrs.PlateCarree(),cmap = cmp1, norm= norm, extend='both')
    
    level_min = 100
    hk = h_isolev[:,:,k[ii]]
    hk[abs(lat)<5] = np.nan
    level_max = np.floor((np.nanmax(hk))/100)*100
    
    if ii < 3:
        levels = np.linspace(level_min,level_max,3)
    else: 
        levels = np.linspace(level_min,level_max,5)
    contour_lines = ax0.contour(lon,lat,np.squeeze(hk),levels,transform = ccrs.PlateCarree(),colors='black',zorder=10000)
    labels = ax0.clabel(
        contour_lines,
        inline=True,
        fontsize=10,
        fmt='%d'
    )
    
    for label in labels:
        label.set_bbox(dict(
            facecolor='white',
            edgecolor='none',
            alpha=0.7,
            boxstyle='round,pad=0.2',
            zorder = 1000
        ))
    
    # Plot MLD mask
    ax0.contourf(
        lon, lat, np.roll(mldcont[:,:,k[ii]],circshift_lim,0),
        transform=ccrs.PlateCarree(),
        colors='black',alpha=0.5,
        levels=[-0.5,0.5],
        zorder=10000,
    )
    
    ax0.contour(
        lon, lat, np.roll(mldcont[:,:,k[ii]],circshift_lim,0),
        transform=ccrs.PlateCarree(),
        colors='black',
        levels=[-0.5,0.5],
        zorder=100000,
    )
    
    land = cfeature.NaturalEarthFeature('physical', 'land', scale='50m', edgecolor='none', facecolor=cfeature.COLORS['land'], linewidth=.25,zorder = 100)
    ax0.add_feature(land, facecolor='k')
    
    # Thicken axis borders
    for spine in ax0.spines.values():
            spine.set_linewidth(2)

# Colorbar 
scale        = 0.1
bounds       = [-10,-8,-6,-4,-2,-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1,2,4,6,8,10]
bounds2      = np.multiply(bounds,scale)
boundticks2  = np.multiply(list(np.arange(-10,12,2)),scale)
norm_cont    = mpl.colors.BoundaryNorm(bounds2, cmp1.N)

cbar = plt.colorbar(
        mpl.cm.ScalarMappable(cmap=cmp1, norm = norm_cont),
    extend='both',
    ticks=boundticks2,
    spacing='proportional',
    shrink=0.5,
    ax=axs,
    location='bottom',
    pad = 0.03,
)

cbar.ax.tick_params(labelsize=18)
cbar.set_label(label=r'$m$ $day^{-1}$',size=20)

plt.savefig('...\figures\figure1.png',bbox_inches='tight', dpi=300)
