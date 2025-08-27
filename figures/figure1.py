# FIGURE 1: Time-mean OLIV3 sigma 26 
import matplotlib.pyplot as plt
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
import cartopy.mpl.ticker as cticker
from cartopy.util import add_cyclic_point
import matplotlib.ticker as mticker
import matplotlib.patches as mpatches
import xarray as xr

anno = np.arange(1993,2016)

filo            = '...\OLIV3\oliv3_glob_025_ekf_annual_isolevm_5filtr.nc'
file            = xr.open_dataset(filo)
lon, lat, iso   = file.variables['longitude'].values.T, file.variables['latitude'].values.T, file.variables['isolev'].values.T
w               = np.mean(file.variables['w_oliv3_isolevf'].values.T,3)

dimlon, dimlat, dimdep = lon.shape[0], lat.shape[1], len(iso)

# Loading Ekman pumping:
    
filo            = '...\ERA5\ekman_pumping_glob_025_filtered_annual.nc'
file            = xr.open_dataset(filo)
lone, late      = file.variables['longitude'].values.T, file.variables['latitude'].values.T
wek             = np.mean(file.variables['wek_era5'].values.T,2)

wek[abs(late)<5] = np.nan

# Colormap
scale           = 10**(-6)
contens_ticks   = np.arange(-10,12,2)*scale
contens         = np.array([-10, -8, -6, -4, -2, -1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1, 2, 4, 6, 8, 10])*scale
norm            = mpl.colors.BoundaryNorm(contens, cmp1.N)

k               = 20 # 20 sigma 26
proj            = ccrs.Robinson(central_longitude=-60)

# Figure:
fig, axs        = plt.subplots(1, 1, figsize=(17, 10),subplot_kw={"projection": ccrs.Robinson(central_longitude=-60)}, layout='constrained',)

ax0 = axs
ax0.set_title(r'Time-mean OLIV3 at $\sigma$ '+str(round(iso[k],2))+ ' kg m$^{-3}$',fontsize = 26)
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

wk = w[:,:,k]
wk[abs(lat)<5] = np.nan
cs = ax0.contourf(lon, lat, wk,contens,transform=ccrs.PlateCarree(),cmap = cmp1, norm= norm,extend='both')

contour_lines = axs.contour(lon,lat,np.squeeze(h_isolev[:,:,k]*mm_eq*maskcc),np.arange(100,600,100),transform = ccrs.PlateCarree(),colors='black',zorder=10000)
labels = axs.clabel(
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
        boxstyle='round,pad=0.2'
    ))

# Plot MLD mask
ax0.contourf(LONimforth, LATimforth, mldcont[:,:,k], transform=ccrs.PlateCarree(), colors='black',alpha=0.5, levels=[-0.5,0.5], zorder=10000)
ax0.contour(LONimforth, LATimforth, mldcont[:,:,k], transform=ccrs.PlateCarree(), colors='black', levels=[-0.5,0.5], zorder=100000)

# Land
land = cfeature.NaturalEarthFeature('physical', 'land', scale='50m', edgecolor='none', facecolor=cfeature.COLORS['land'], linewidth=.25,zorder = 100)
ax0.add_feature(land, facecolor='k')

# Colorbar ----------------------------------------------------------------
cbar = fig.colorbar(mpl.cm.ScalarMappable(cmap=cmp1, norm=norm2), extend='both', ticks=boundticks, spacing='proportional', orientation='vertical', shrink=0.7, ax = ax0)

cbar.ax.tick_params(labelsize=20)
cbar.set_label(label='Vertical velocities [10$^{-6}$ m/s]',size = 20)
cbar.ax.set_yticklabels(boundticks)

# Thicken axis borders
for spine in ax0.spines.values():
        spine.set_linewidth(2)

# Rectangles Figure 2:
ax0.add_patch(mpatches.Rectangle(xy=[-30, 8], width=14, height=8,facecolor='none',edgecolor='darkviolet',linewidth=2,transform=ccrs.PlateCarree(),zorder = 10000))
ax0.add_patch(mpatches.Rectangle(xy=[-50, 15], width=20, height=15,facecolor='none',edgecolor='darkviolet',linewidth=2,transform=ccrs.PlateCarree(),zorder = 10000))
ax0.add_patch(mpatches.Rectangle(xy=[-55, 50], width=25, height=10,facecolor='none',edgecolor='darkviolet',linewidth=2,transform=ccrs.PlateCarree(),zorder = 10000))
        
plt.savefig('...\figures\figure1_VF.png',bbox_inches='tight', dpi=300)



