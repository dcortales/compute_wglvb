# Figure 2: Hovmoller diagram OLIV3

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

filo            = '...\OLIV3\oliv3_glob_025_ekf_annual_isolevm_5filtr.nc'
file            = xr.open_dataset(filo)

filok           = '...\ERA5\ekman_pumping_glob_025_filtered_annual.nc'
filek           = xr.open_dataset(filok)
lone, late      = filek.variables['longitude'].values.T, filek.variables['latitude'].values.T
wek             = filek.variables['wek_era5'].values.T

wek[abs(late)<5] = np.nan

# Regionalisation:

#TG:
min_lat, max_lat = 8, 16
min_lon, max_lon = 360-30, 360-16
min_iso, max_iso = 18,29
    
trop_mask                                       = np.ones(lon.shape)
trop_mask[(lon < min_lon) | (lon > max_lon)]    = np.nan
trop_mask[(lat < min_lat) | (lat > max_lat)]    = np.nan
wt = np.zeros(w.shape[2])
for k in np.arange(0,w.shape[2]):
    wt[k]= np.nanmean(w[:,:,k]*trop_mask)
    
wekt = np.zeros(len(anno))
for t in np.arange(0,len(anno)):
    wekt[t] = np.nanmean(wek[:,:,t]*trop_mask)
    

w_isopt = np.zeros((len(iso),len(anno)))
for t in np.arange(0,len(anno)):
    w   = file['w_oliv3_isolevf'].isel(time=t).values.T
    for k in np.arange(0,len(iso)):
        if k == 0:
            wk = w[:,:,k]*trop_mask
            wkk = wk[np.isnan(wk) == False]
            w_isopt[k,t] = np.nanmean(wkk)
        else:
            wk = w[:,:,k]*trop_mask
            wkk = wk[np.isnan(wk) == False]
            w_isopt[k,t] = np.nanmean(wkk)

w_isopt[13,:] = wekt
    
#STG:
min_lat, max_lat = 15, 30
min_lon, max_lon = 360-65, 360-45
min_iso, max_iso = 18,29
    
subt_mask = np.ones(lon.shape)
subt_mask[(lon < min_lon)| (lon > max_lon)] = np.nan
subt_mask[(lat < min_lat) | (lat > max_lat)] = np.nan

wt = np.zeros(w.shape[2])
for k in np.arange(0,w.shape[2]):
    wt[k]= np.nanmean(w[:,:,k]*subt_mask)
    
weks = np.zeros(len(anno))
for t in np.arange(0,len(anno)):
    weks[t] = np.nanmean(wek[:,:,t]*subt_mask)
    
w_isops = np.zeros((len(iso),len(anno)))
for t in np.arange(0,len(anno)):
    w   = file['w_oliv3_isolevf'].isel(time=t).values.T
    for k in np.arange(0,len(iso)):
        if k == 0:
            wk = w[:,:,k]*subt_mask
            wkk = wk[np.isnan(wk) == False]
            w_isops[k,t] = np.nanmean(wkk)
        else:
            wk = w[:,:,k]*subt_mask
            wkk = wk[np.isnan(wk) == False]
            w_isops[k,t] = np.nanmean(wkk)

w_isops[18,:] = weks

#SPG:
min_lat, max_lat = 50,60
min_lon, max_lon = 360-55, 360-30
min_iso, max_iso = 18,29
    
subp_mask = np.ones(lon.shape)
subp_mask[(lon < min_lon)| (lon > max_lon)] = np.nan
subp_mask[(lat < min_lat) | (lat > max_lat)] = np.nan

wt = np.zeros(w.shape[2])
for k in np.arange(0,w.shape[2]):
    wt[k]= np.nanmean(w[:,:,k]*subp_mask)

wekp = np.zeros(len(anno))
for t in np.arange(0,len(anno)):
    wekp[t] = np.nanmean(wek[:,:,t]*subp_mask)

w_isopp = np.zeros((len(iso),len(anno)))
for t in np.arange(0,len(anno)):
    w   = file['w_oliv3_isolevf'].isel(time=t).values.T
    for k in np.arange(0,len(iso)):
        if k == 0:
            wk = w[:,:,k]*subp_mask
            wkk = wk[np.isnan(wk) == False]
            w_isopp[k,t] = np.nanmean(wkk)
        else:
            wk = w[:,:,k]*subp_mask
            wkk = wk[np.isnan(wk) == False]
            w_isopp[k,t] = np.nanmean(wkk)

w_isopp[35,:] = wekp

# Figure:

scale       = 10**(-6)
vmin = -1.5e-6
vmax = 1.5e-6
cmp = 'RdBu_r'

fig, axs = plt.subplots(3, 1, figsize=(7, 12))

# Panel 1: TG
pc0 = axs[0].pcolor(anno, iso[14:22], w_isopt[14:22,:], cmap=cmp, vmin=vmin, vmax=vmax)
contours = axs[0].contour(anno, iso[13:23], w_isopt[13:23,:]/scale, colors='black', vmin=vmin/scale, vmax=vmax/scale)
plt.clabel(contours, inline=True, fontsize=8)
axs[0].invert_yaxis()
axs[0].set_title('(a) NATL Tropical gyre [30-16$^o$W]:[8-16$^o$N]',fontsize=16)
axs[0].set_ylabel('Isopycnal level (kg/m$^3$)',fontsize=14)
axs[0].set_xticklabels([])
axs[0].set_xlim([1992.5,2015.5])
axs[0].set_ylim([(iso[21]+iso[22])/2,(iso[13]+iso[14])/2])

# Panel 2: STG
pc1 = axs[1].pcolormesh(anno, iso[19:29], w_isops[19:29,:], cmap=cmp, vmin=vmin, vmax=vmax)
contours = axs[1].contour(anno, iso[18:30], w_isops[18:30,:]/scale, colors='black', vmin=vmin/scale, vmax=vmax/scale)
plt.clabel(contours, inline=True, fontsize=8)
axs[1].invert_yaxis()
axs[1].set_title('(b) NATL Subtropical gyre [50-30$^o$W]:[15-30$^o$N]',fontsize=16)
axs[1].set_ylabel('Isopycnal level (kg/m$^3$)',fontsize=14)
axs[1].set_xticklabels([])
axs[1].set_xlim([1992.5,2015.5])
axs[1].set_ylim([(iso[28]+iso[29])/2,(iso[18]+iso[19])/2])

# Panel 3: SPG
pc2 = axs[2].pcolormesh(anno, iso[36:55], w_isopp[36:55,:], cmap=cmp, vmin=vmin, vmax=vmax)
contours = plt.contour(anno, iso[35:56], w_isopp[35:56,:]/scale, colors='black', vmin=vmin/scale, vmax=vmax/scale)
plt.clabel(contours, inline=True, fontsize=8)
axs[2].invert_yaxis()
axs[2].set_title('(c) NATL Subpolar gyre [55-30$^o$W]:[50-60$^o$N]',fontsize=16)
axs[2].set_ylabel('Isopycnal level (kg/m$^3$)',fontsize=14)
axs[2].set_xlabel('Years',fontsize=14)
axs[2].set_xlim([1992.5,2015.5])
axs[2].set_ylim([(iso[54]+iso[55])/2,(iso[35]+iso[36])/2])

for ax in axs:
    for spine in ax.spines.values():
        spine.set_linewidth(2)
    ax.tick_params(axis='both', width=2, length=6, labelsize=12)
        
pos = axs[2].get_position()
cax = fig.add_axes([pos.x0, 0.95, pos.width, 0.015])
scale       = 10**(-7)
contens_ticks     = np.arange(-15,17.5,2.5)*scale
boundticks  = np.round(np.arange(-1.5,1.75,0.25),2)
cb = fig.colorbar(pc2,cax=cax, orientation='horizontal',extend='both',ticks=contens_ticks)
cb.set_label('Vertical velocity [10$^{-6}$ m/s)', labelpad=-60, fontsize=14, loc='center')
cb.ax.set_xticklabels(boundticks,fontsize=12)

plt.savefig('...\figures\figure2_VF.png',bbox_inches='tight', dpi=300)
