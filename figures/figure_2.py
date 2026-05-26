# FIGURE 1: Time-mean OLIV3 sigma 26 and Ekman pumping
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
import scipy.io as sio
from matplotlib.colors import ListedColormap
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.patches as mpatches
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

# Ekman pumping
filo            = '...\Data\ERA5\ekman_pumping_glob_025_filtered_annual.nc'
file            = xr.open_dataset(filo)
lone, late      = file.variables['longitude'].values.T, file.variables['latitude'].values.T
wek             = np.mean(file.variables['wek_era5'].values.T,2)

wek[abs(late)<5] = np.nan

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
k           = 20 # 20 sigma 26

fig, axs    = plt.subplots(2, 1, figsize=(17, 17),subplot_kw={"projection": ccrs.Robinson(central_longitude=-60)}, layout='constrained',)

# Panel (a): OLIV3 sigma 26
ax0 = axs[0]
ax0.set_title(r'(a) Time-mean OLIV3 at $\sigma$ '+str(round(iso[k],2))+ ' kg m$^{-3}$',fontsize = 24)
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
grd.xlabel_style = {'size': 20}
grd.ylabel_style = {'size': 20}  

# Plot oliv3
wk = w[:,:,k]
wk[abs(lat)<5] = np.nan
cs = ax0.contourf(lon, lat, wk*60*60*24,bounds2,transform=ccrs.PlateCarree(),cmap = cmp1, norm= norm, extend='both')

# Mask h_lines
maskcc                      = np.zeros(wk.shape)
maskcc[np.isnan(wk)==True]  = np.nan
maskcc[np.isnan(wk)==False] = 1
mm_eq                       = np.ones(maskcc.shape)
mm_eq[abs(lat)<5]           = np.nan

# Plot depth contour lines
contour_lines   = ax0.contour(lon,lat,np.squeeze(h_isolev[:,:,k]*mm_eq*maskcc),np.arange(100,600,100),transform = ccrs.PlateCarree(),colors='black',zorder=10000)
labels          = ax0.clabel(contour_lines, inline=True, fontsize=10, fmt='%d')

for label in labels:
    label.set_bbox(dict(facecolor='white', edgecolor='none', alpha=0.7, boxstyle='round,pad=0.2'))

# Plot MLD mask
ax0.contourf(lon, lat, np.roll(mldcont[:,:,k],circshift_lim,0), transform=ccrs.PlateCarree(), colors='black',alpha=0.5, levels=[-0.5,0.5], zorder=10000,)
ax0.contour(lon, lat, np.roll(mldcont[:,:,k],circshift_lim,0), transform=ccrs.PlateCarree(), colors='black', levels=[-0.5,0.5], zorder=100000,)

# Land
land = cfeature.NaturalEarthFeature('physical', 'land', scale='50m', edgecolor='none', facecolor=cfeature.COLORS['land'], linewidth=.25,zorder = 100)
ax0.add_feature(land, facecolor='k')

for spine in ax0.spines.values():
        spine.set_linewidth(2)
        
# Rectangles Fig. 3:
ax0.add_patch(mpatches.Rectangle(xy=[-30, 8], width=14, height=8,
                                    facecolor='none',
                                    edgecolor='darkviolet',
                                    linewidth=2,
                                    transform=ccrs.PlateCarree(),
                                    zorder = 10000))

ax0.add_patch(mpatches.Rectangle(xy=[-50, 15], width=20, height=15,
                                    facecolor='none',
                                    edgecolor='darkviolet',
                                    linewidth=2,
                                    transform=ccrs.PlateCarree(),
                                    zorder = 10000))

ax0.add_patch(mpatches.Rectangle(xy=[-55, 50], width=25, height=10,
                                    facecolor='none',
                                    edgecolor='darkviolet',
                                    linewidth=2,
                                    transform=ccrs.PlateCarree(),
                                    zorder = 10000))

# Panel (b): Ekman pumping
ax1 = axs[1]
ax1.set_title(r'(b) Time-mean ERA5 Ekman pumping',fontsize = 24)
ax1.set_xlabel(r'$Longitude(^o)$')
ax1.set_ylabel(r'$Latitude(^o)$')

grd = ax1.gridlines(
    draw_labels=True, 
    xlocs=range(-180, 181, 90), 
    ylocs=range(-60, 61, 30), 
    color='gray',
    zorder=10000,
)
grd.top_labels = False
grd.xlabel_style = {'size': 18}
grd.ylabel_style = {'size': 18}

cs = ax1.contourf(lon, lat, wek*60*60*24,bounds2,transform=ccrs.PlateCarree(),cmap = cmp1, norm= norm)

# sigma 26 mask
mask26 = np.ones(wk.shape)
mask26[np.isnan(wk)==False] = np.nan
mask26[np.isnan(wek)==True] = np.nan

plt.rcParams['hatch.color'] = (1,1,1)
plt.rcParams['hatch.linewidth'] = 3
ax1.contourf(lon, lat, mask26,[0,1],transform = ccrs.PlateCarree(),hatches=['///'],colors='none',extend='both',zorder=1000)

land = cfeature.NaturalEarthFeature('physical', 'land', scale='50m', edgecolor='none', facecolor=cfeature.COLORS['land'], linewidth=.25,zorder = 100)
ax1.add_feature(land, facecolor='k')

scale        = 0.1
bounds       = [-10,-8,-6,-4,-2,-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1,2,4,6,8,10]
bounds2      = np.multiply(bounds,scale)
boundticks2  = np.multiply(list(np.arange(-10,12,2)),scale)
norm_cont    = mpl.colors.BoundaryNorm(bounds2, cmp1.N)

# Colorbar
cb = fig.colorbar(
    mpl.cm.ScalarMappable(cmap=cmp1, norm=norm),
    extend='both',
    ticks=boundticks2,
    spacing='proportional',
    orientation='vertical',
    label='Relative error [%]',
    shrink=0.4,
    ax=axs.ravel().tolist()
)

cb.ax.tick_params(labelsize=18)
cb.set_label(label='$m$ $day^{-1}$',size = 18)

for spine in ax1.spines.values():
        spine.set_linewidth(2)

plt.savefig('...\figures\figure2.png',bbox_inches='tight', dpi=300)
