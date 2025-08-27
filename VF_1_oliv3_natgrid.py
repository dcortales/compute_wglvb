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


# GRID loading:
latwm           = sio.loadmat('C:\\Users\\yago_\\Documents\\LOCEAN\\OLIV3 paper\\LVBw3D_ARMOR3D_GLOB.mat')
LONimforth      = latwm['LONimone']

LATimforth      = latwm['LATimone']
isoplm          = sio.loadmat('w_isop_distn.mat')
isopl           = isoplm['isopl']
lonwm           = sio.loadmat('C:\\Users\\yago_\\Documents\\LOCEAN\\OLIV3 paper\\Limone_GLOB.mat')
LONimone        = lonwm['LONimone']
latwm           = sio.loadmat('C:\\Users\\yago_\\Documents\\LOCEAN\\OLIV3 paper\\Limone_GLOB.mat')
LATimone        = latwm['LATimone']

# Variables loading: 
    
wofintm         = sio.loadmat('C:\\Users\\yago_\\Documents\\LOCEAN\\OLIV3 paper\\wofint_masked.mat')
wof             = wofintm['wofint_masked']

MLDm            = sio.loadmat('C:\\Users\\yago_\\Documents\\LOCEAN\\OLIV3 paper\\MLD_ARMOR3D_GLOB.mat')
mldcont         = MLDm['mld_contac']

# LVB error:
eem         = sio.loadmat('C:\\Users\\yago_\\Documents\\LOCEAN\\OLIV3 paper\\LVB_ee_GLOB_int.mat')
eelvb       = eem['eeint']

# Colormap: 
    
RBensm      = sio.loadmat('figure1/RBens.mat')
RBens       = RBensm['RBens']
RBens1m     = sio.loadmat('figure1/RBens1.mat')
RBens1      = RBens1m['RBens1']

cmp         = ListedColormap(RBens)
cmp1        = ListedColormap(RBens1)

# Circshifted

lonwm           = sio.loadmat('C:\\Users\\yago_\\Documents\\LOCEAN\\OLIV3 paper\\ww_ARMOR3D_isopl_filtr_distn_masked.mat')
LONim           = lonwm['LONfm']
latwm           = sio.loadmat('C:\\Users\\yago_\\Documents\\LOCEAN\\OLIV3 paper\\ww_ARMOR3D_isopl_filtr_distn_masked.mat')
LATim           = latwm['LATm']
wfcomp_isopm    = sio.loadmat('C:\\Users\\yago_\\Documents\\LOCEAN\\OLIV3 paper\\ww_ARMOR3D_isopl_filtr_distn_masked.mat')
w_isop          = wfcomp_isopm['ww_isop_distnm']


# Bounds:

scale       = 10**(-6)
contens     = [-200, -10, -8, -6, -4, -2, -1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1, 2, 4, 6, 8, 10, 200]
contens2    = np.multiply(contens,scale)
cont4       = np.array([-20, -1, -0.8, -0.6, -0.2, -0.05, 0, 0.05, 0.2, 0.6, 0.8, 1])*100;
ca          = [-4000, -40, -30, -20, -10, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 10, 20, 30, 40, 400]

bounds      = [-1,-0.8,-0.6,-0.4,-0.2,-0.1,-0.08,-0.06,-0.04,-0.02,0,0.02,0.04,0.06,0.08,0.1,0.2,0.4,0.6,0.8,1]
bounds2     = [-10,-8,-6,-4,-2,-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1,2,4,6,8,10]
boundticks  = [-10,-8,-6,-4,-2,0,2,4,6,8,10]
norm        = mpl.colors.BoundaryNorm(bounds, cmp1.N)
norm2       = mpl.colors.BoundaryNorm(bounds2, cmp1.N)

k = 20

# LVB_mask ------------------------------------------------------------------
LVB_mask                            = np.zeros_like(eelvb[:,:,k])
LVB_mask[LVB_mask==0]               = np.nan
LVB_mask[abs(eelvb[:,:,k])>10]      = 1
LVB_mask[abs(LATimone)<5]           = np.nan

# obs. w_{GLVB} figure -----------------------------------------------------

proj = ccrs.Robinson(central_longitude=-60)

fig = plt.figure(figsize=(20, 10))
ax0 = plt.axes(projection=proj)
ax0.set_title(r'$OLIV3$ at $\sigma$ '+str(round(isopl[0,k],2))+ ' kg m$^{-3}$',fontsize = 22)
ax0.set_xlabel(r'$Longitude(^o)$')
ax0.set_ylabel(r'$Latitude(^o)$')

grd = ax0.gridlines(
    draw_labels=True, 
    xlocs=range(-180, 181, 90), 
    ylocs=range(-60, 61, 30), 
    color='gray',
    zorder=10000,
)
grd.top_labels      = False
grd.xlabel_style    = {'size': 18} 
grd.ylabel_style    = {'size': 18}    

w_isopk = w_isop[:,:,k]
w_isopk.transpose()
ax0.contourf(LONim, LATim, w_isopk,contens2,transform = ccrs.PlateCarree(),cmap = cmp, vmin = -10*10**(-6), vmax = 10*10**(-6),zorder=10)
#ax0.contour(LONimone, LATimone, wof[:,:,k],[-1,0,1],transform = ccrs.PlateCarree(),colors = 'k', vmin = -1, vmax = 1,zorder=100)

# Colorbar ----------------------------------------------------------------
cbar = fig.colorbar(
        mpl.cm.ScalarMappable(cmap=cmp1, norm=norm2),
    extend='both',
    ticks=boundticks,
    spacing='proportional',
    orientation='vertical',
    shrink=0.8,
    ax = ax0
)

cbar.ax.tick_params(labelsize=18)
cbar.set_label(label=r'(10$^{-6}$) m s$^{-1}$',size=20)   

# LVB mask ----------------------------------------------------------------
plt.rcParams['hatch.color'] = (0,0,0)
plt.rcParams['hatch.linewidth'] = 1.2
#ax0.contourf(LONimone, LATimone, LVB_mask,[0,1],transform = ccrs.PlateCarree(),hatches=['///'],colors='none',extend='both',zorder=1000)


land = cfeature.NaturalEarthFeature('physical', 'land', scale='50m', edgecolor='none', facecolor=cfeature.COLORS['land'], linewidth=.25,zorder = 100)
ax0.add_feature(land, facecolor='k')

# %% Correlation coefficient ERA5 and OLIV3

filo            = 'C:\\Users\yago_\Documents\LOCEAN\Data\OLIV3\oliv3_glob_025_ekf_annual_isolevm_5filtr.nc'
file            = xr.open_dataset(filo)
lon, lat, iso   = file.variables['longitude'].values.T, file.variables['latitude'].values.T, file.variables['isolev'].values.T
wt              = file.variables['w_oliv3_isolevf'].values.T

dimlon, dimlat, dimdep = lon.shape[0], lat.shape[1], len(iso)

# Loading Ekman pumping:
    
filo            = 'C:\\Users\yago_\Documents\LOCEAN\Data\ERA5\ekman_pumping_glob_025_filtered_annual.nc'
file            = xr.open_dataset(filo)
lone, late      = file.variables['longitude'].values.T, file.variables['latitude'].values.T
wekt            = file.variables['wek_era5'].values.T

k = 20
cc = np.zeros(lon.shape)
for i in np.arange(0,cc.shape[0]):
    for j in np.arange(0,cc.shape[1]):
        R = np.corrcoef(wekt[i,j,:],wt[i,j,k,:])
        cc[i,j] = R[0,1]

# %% Correlation figure

filo            = 'C:\\Users\yago_\Documents\LOCEAN\Data\OLIV3\oliv3_glob_025_ekf_annual_isolevm.nc'
file            = xr.open_dataset(filo)
h_isolev        = file.variables['h_isolev'].values.T
h_isolev = h_isolev[:,:,:,0]
# %%
proj        = ccrs.Robinson(central_longitude=-60)

fig, axs    = plt.subplots(1, 1, figsize=(17, 10),subplot_kw={"projection": ccrs.Robinson(central_longitude=-60)}, layout='constrained',)

BrBGcentre  = plt.get_cmap('RdYlBu_r')
cont_rho    = [-1, 0, 0.5,0.6,0.7,0.8,0.9,1]
boundstick  = np.arange(-1,1.2,0.2)
bounds      = [-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]
norm        = mpl.colors.BoundaryNorm(bounds, BrBGcentre.N)

k = 20

#plt.figure(figsize=(20, 10))
axs.set_title(r'Correlation coefficient between ERA5 $w_{Ek}$ and OLIV3 at $\sigma$ '+str(round(isopl[0,k],2))+ ' kg m$^{-3}$',fontsize = 24)
axs.set_xlabel(r'$Longitude(^o)$')
axs.set_ylabel(r'$Latitude(^o)$')

grd = axs.gridlines(
    draw_labels=True, 
    xlocs=range(-180, 181, 90), 
    ylocs=range(-60, 61, 30), 
    color='gray',
    zorder=10000,
)
grd.top_labels = False
grd.xlabel_style = {'size': 18}
grd.ylabel_style = {'size': 18}

axs.contourf(lon, lat, cc,bounds,transform = ccrs.PlateCarree(),cmap = BrBGcentre, norm = norm,zorder=10)

maskcc = np.zeros(cc.shape)
maskcc[np.isnan(cc)==True] = np.nan
maskcc[np.isnan(cc)==False] = 1

mm_eq = np.ones(maskcc.shape)
mm_eq[abs(lat)<5] = np.nan

contour_lines = axs.contour(lon,lat,np.squeeze(h_isolev[:,:,k]*mm_eq*maskcc),np.arange(100,600,100),transform = ccrs.PlateCarree(),colors='black',zorder=100)

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

land = cfeature.NaturalEarthFeature('physical', 'land', scale='50m', edgecolor='none', facecolor=cfeature.COLORS['land'], linewidth=.25,zorder = 100)
axs.add_feature(land, facecolor='k')

#plt.figure(figsize=(5, 5))
#ax0 = plt.axes(projection=ccrs.PlateCarree())

cbar = plt.colorbar(
        mpl.cm.ScalarMappable(cmap=BrBGcentre,norm=norm),
    extend='both',
    ticks=boundstick,
    spacing='proportional',
    orientation='vertical',
    shrink=0.8,
    ax = axs
)

cbar.ax.tick_params(labelsize=17)
cbar.set_label(label='Correlation coefficient',size = 18)

for spine in axs.spines.values():
    spine.set_linewidth(2)

plt.savefig('1_oliv3_natgrid_corrcoeff_VF.png',bbox_inches='tight', dpi=300)

# %% OLIV3 sigma 26

import numpy as np
import xarray as xr
from netCDF4 import Dataset as datas
import h5py

import matplotlib.pyplot as plt

anno = np.arange(1993,2016)

filo            = 'C:\\Users\yago_\Documents\LOCEAN\Data\OLIV3\oliv3_glob_025_ekf_annual_isolevm_5filtr.nc'
file            = xr.open_dataset(filo)
lon, lat, iso   = file.variables['longitude'].values.T, file.variables['latitude'].values.T, file.variables['isolev'].values.T
w               = np.mean(file.variables['w_oliv3_isolevf'].values.T,3)

dimlon, dimlat, dimdep = lon.shape[0], lat.shape[1], len(iso)

# Loading Ekman pumping:
    
filo            = 'C:\\Users\yago_\Documents\LOCEAN\Data\ERA5\ekman_pumping_glob_025_filtered_annual.nc'
file            = xr.open_dataset(filo)
lone, late      = file.variables['longitude'].values.T, file.variables['latitude'].values.T
wek             = np.mean(file.variables['wek_era5'].values.T,2)

wek[abs(late)<5] = np.nan

# Colormap

scale       = 10**(-6)
contens_ticks     = np.arange(-10,12,2)*scale
contens     = np.array([-10, -8, -6, -4, -2, -1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1, 2, 4, 6, 8, 10])*scale
norm        = mpl.colors.BoundaryNorm(contens, cmp1.N)

k           = 20 # 20 sigma 26
proj        = ccrs.Robinson(central_longitude=-60)

# %%
fig, axs    = plt.subplots(1, 1, figsize=(17, 10),subplot_kw={"projection": ccrs.Robinson(central_longitude=-60)}, layout='constrained',)

# oliv3:
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
ax0.contourf(
    LONimforth, LATimforth, mldcont[:,:,k],
    transform=ccrs.PlateCarree(),
    colors='black',alpha=0.5,
    levels=[-0.5,0.5],
    zorder=10000,
)

ax0.contour(
    LONimforth, LATimforth, mldcont[:,:,k],
    transform=ccrs.PlateCarree(),
    colors='black',
    levels=[-0.5,0.5],
    zorder=100000,
)

land = cfeature.NaturalEarthFeature('physical', 'land', scale='50m', edgecolor='none', facecolor=cfeature.COLORS['land'], linewidth=.25,zorder = 100)
ax0.add_feature(land, facecolor='k')

# Colorbar ----------------------------------------------------------------
cbar = fig.colorbar(
        mpl.cm.ScalarMappable(cmap=cmp1, norm=norm2),
    extend='both',
    ticks=boundticks,
    spacing='proportional',
    orientation='vertical',
    shrink=0.7,
    ax = ax0
)

cbar.ax.tick_params(labelsize=20)
cbar.set_label(label='Vertical velocities [10$^{-6}$ m/s]',size = 20)
cbar.ax.set_yticklabels(boundticks)

# Thicken axis borders
for spine in ax0.spines.values():
        spine.set_linewidth(2)
        

ax0.add_patch(mpatches.Rectangle(xy=[-30, 8], width=14, height=8,
                                    facecolor='none',
                                    #alpha=0.2,
                                    edgecolor='darkviolet',
                                    linewidth=2,
                                    transform=ccrs.PlateCarree(),
                                    zorder = 10000))

ax0.add_patch(mpatches.Rectangle(xy=[-50, 15], width=20, height=15,
                                    facecolor='none',
                                    #alpha=0.2,
                                    edgecolor='darkviolet',
                                    linewidth=2,
                                    transform=ccrs.PlateCarree(),
                                    zorder = 10000))

ax0.add_patch(mpatches.Rectangle(xy=[-55, 50], width=25, height=10,
                                    facecolor='none',
                                    #alpha=0.2,
                                    edgecolor='darkviolet',
                                    linewidth=2,
                                    transform=ccrs.PlateCarree(),
                                    zorder = 10000))
        
plt.savefig('1_oliv3_natgrid_VF.png',bbox_inches='tight', dpi=300)



# %% Hovmoller diagram OLIV3 + Ekman pumping ---------------------------------

filo            = 'C:\\Users\yago_\Documents\LOCEAN\Data\OLIV3\oliv3_glob_025_ekf_annual_isolevm_5filtr.nc'
file            = xr.open_dataset(filo)

filok            = 'C:\\Users\yago_\Documents\LOCEAN\Data\ERA5\ekman_pumping_glob_025_filtered_annual.nc'
filek            = xr.open_dataset(filok)
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

## Correlation coeficient

Rt = np.zeros(len(iso))
for k in np.arange(0,len(iso)):
    R = np.corrcoef(wekt,w_isopt[k,:])
    Rt[k] = R[0,1]

Rs = np.zeros(len(iso))
for k in np.arange(0,len(iso)):
    R = np.corrcoef(weks,w_isops[k,:])
    Rs[k] = R[0,1]

Rp = np.zeros(len(iso))
for k in np.arange(0,len(iso)):
    R = np.corrcoef(wekp,w_isopp[k,:])
    Rp[k] = R[0,1] 
    
## RMSE

RMSEt, RMSEs, RMSEp = np.zeros(len(iso)), np.zeros(len(iso)), np.zeros(len(iso))
for k in np.arange(0,len(iso)):
    RMSEt[k] = np.sqrt(np.mean((wekt-w_isopt[k,:])**2))
    RMSEs[k] = np.sqrt(np.mean((weks-w_isops[k,:])**2))
    RMSEp[k] = np.sqrt(np.mean((wekp-w_isopp[k,:])**2))

## FIGURE

scale       = 10**(-6)
 
vmin = -1.5e-6
vmax = 1.5e-6
cmp = 'RdBu_r'

# Create subplots
fig, axs = plt.subplots(3, 1, figsize=(7, 12))

# Panel 1: Original variable with pcolormesh
pc0 = axs[0].pcolor(anno, iso[14:22], w_isopt[14:22,:], cmap=cmp, vmin=vmin, vmax=vmax)
contours = axs[0].contour(anno, iso[13:23], w_isopt[13:23,:]/scale, colors='black', vmin=vmin/scale, vmax=vmax/scale)
plt.clabel(contours, inline=True, fontsize=8)
axs[0].invert_yaxis()
axs[0].set_title('(a) NATL Tropical gyre [30-16$^o$W]:[8-16$^o$N]',fontsize=16)
axs[0].set_ylabel('Isopycnal level (kg/m$^3$)',fontsize=14)
axs[0].set_xticklabels([])
#axs[0,0].hlines((iso[13]+iso[14])/2,1992,2020,color='black')
axs[0].set_xlim([1992.5,2015.5])
axs[0].set_ylim([(iso[21]+iso[22])/2,(iso[13]+iso[14])/2])


# Panel 2: First additional variable
pc1 = axs[1].pcolormesh(anno, iso[19:29], w_isops[19:29,:], cmap=cmp, vmin=vmin, vmax=vmax)
contours = axs[1].contour(anno, iso[18:30], w_isops[18:30,:]/scale, colors='black', vmin=vmin/scale, vmax=vmax/scale)
plt.clabel(contours, inline=True, fontsize=8)
axs[1].invert_yaxis()
axs[1].set_title('(b) NATL Subtropical gyre [50-30$^o$W]:[15-30$^o$N]',fontsize=16)
axs[1].set_ylabel('Isopycnal level (kg/m$^3$)',fontsize=14)
axs[1].set_xticklabels([])
#axs[1,0].hlines((iso[18]+iso[19])/2,1992,2020,color='black')
axs[1].set_xlim([1992.5,2015.5])
axs[1].set_ylim([(iso[28]+iso[29])/2,(iso[18]+iso[19])/2])

# Panel 3: Second additional variable
pc2 = axs[2].pcolormesh(anno, iso[36:55], w_isopp[36:55,:], cmap=cmp, vmin=vmin, vmax=vmax)
contours = plt.contour(anno, iso[35:56], w_isopp[35:56,:]/scale, colors='black', vmin=vmin/scale, vmax=vmax/scale)
plt.clabel(contours, inline=True, fontsize=8)
axs[2].invert_yaxis()
axs[2].set_title('(c) NATL Subpolar gyre [55-30$^o$W]:[50-60$^o$N]',fontsize=16)
axs[2].set_ylabel('Isopycnal level (kg/m$^3$)',fontsize=14)
axs[2].set_xlabel('Years',fontsize=14)
#axs[2,0].hlines((iso[35]+iso[36])/2,1992,2020,color='black')
#fig.colorbar(pc2, ax=axs[2,0], orientation='vertical', label='w [m/s]')
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

plt.savefig('101_oliv3_hovmoller_diag_VF.png',bbox_inches='tight', dpi=300)

