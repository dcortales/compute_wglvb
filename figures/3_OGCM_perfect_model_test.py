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
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.colors as mcolors

file_rho    = sio.loadmat('.../OCCITENS/rho_w_wglvb_occitens_glob_025_annual_isolevm_filtr.mat')
lon         = file_rho['LONimone']
lat         = file_rho['LATimone']
isopl       = file_rho['isopl']
rho         = file_rho['rho_int']

file_relerr = sio.loadmat('.../OCCITENS/relerr_w_wglvb_occitens_glob_025_annual_isolevm_filtr.mat')
relerr      = file_relerr['relerr_int']

file_lvb    = sio.loadmat('.../OCCITENS/lvb_validity_occitens_glob_025_23yr_isolevm_filtr.mat')
lvb         = file_lvb['ee_isop_int']

MLDm        = sio.loadmat('.../OCCITENS/MLDcont_int_GLOB.mat')
mldcont     = MLDm['mld_contint']

# %% Bounds:
bounds      = np.arange(-1,1.2,0.2)
boundticks  = [-1,-0.8,-0.6,-0.4,-0.2,-0.1,0,0.1,0.2,0.4,0.6,0.8,1]

colors      = ['#f72585', '#7209b7']
colors      = ['gold','darkorange', '#47026c']
cmap1       = ListedColormap(colors)

cmap = mpl.cm.gist_earth_r

# Definition truncate colormap:
def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = mcolors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap

# Colormap LVB:
BrBGcentre          = plt.get_cmap('RdYlBu',13)
newcc               = BrBGcentre(np.linspace(0, 1, 13))
yl                  = newcc[6,:]

cmp_rel             = plt.get_cmap('magma_r')
new_cmap_rel        = truncate_colormap(cmp_rel, 0.1, 0.8)

newcolors           = new_cmap_rel(np.linspace(0, 1, 15))
newcolors1          = newcolors
newcolors1[0,:]     = [1, 1, 0.6, 1]
newcmp_rel          = ListedColormap(newcolors1)

# Bounds:
cont_err_rel        = [0,10,30,50,70,90,110,200000]
bounds_rel          = [0,10,30,50,70,90,110]
boundticks_rel      = [0,10,30,50,70,90,110]
norm_rel            = mpl.colors.BoundaryNorm(bounds_rel, newcmp_rel.N)

# Figure 3: OGCM perfect model test

k = 20 # sigma 26
proj = ccrs.Robinson(central_longitude=-60)

fig, axs = plt.subplots(3, 1, figsize=(17, 23),subplot_kw={"projection": ccrs.Robinson(central_longitude=-60)}, layout='constrained',)

ax0 = axs[0]
ax0.set_title(r'(a) Absolute relative error between $w_{g}$ and $w_{tot}$ at $\sigma$ '+str(round(isopl[0,k],2))+ ' kg m$^{-3}$',fontsize = 28)
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
grd.xlabel_style = {'size': 22}
grd.ylabel_style = {'size': 22}  

# Panel a: relative error
cs = ax0.contourf(lon, lat, abs(relerr[:,:,k]),transform=ccrs.PlateCarree(),cmap = newcmp_rel,levels=cont_err_rel,vmin = -0, vmax = 110,zorder=1000)

# LVB mask
lvb_mask = np.zeros(lvb.shape)
lvb_mask = lvb_mask*np.nan
lvb_mask[abs(lvb)>10] = 1
lvb_mask[abs(lat)<4.3] = np.nan
plt.rcParams['hatch.color'] = (0,0,0)
plt.rcParams['hatch.linewidth'] = 1.2
ax0.contourf(lon, lat, lvb_mask[:,:,k],[0,1],transform = ccrs.PlateCarree(),hatches=['///'],colors='none',extend='both',zorder=1000)

# Plot MLD mask
ax0.contourf(lon, lat, mldcont[:,:,k],transform=ccrs.PlateCarree(),colors='black',alpha=0.5,levels=[-0.5,0.5],zorder=100000)
ax0.contour(lon, lat, mldcont[:,:,k],transform=ccrs.PlateCarree(),colors='black',levels=[-0.5,0.5],zorder=100000)

land = cfeature.NaturalEarthFeature('physical', 'land', scale='50m', edgecolor='none', facecolor=cfeature.COLORS['land'], linewidth=.25,zorder = 100)
ax0.add_feature(land, facecolor='k')

cb = fig.colorbar(mpl.cm.ScalarMappable(cmap=newcmp_rel, norm=norm_rel),extend='max',ticks=bounds_rel,spacing='proportional',orientation='vertical',label='Relative error [%]',shrink=0.8,ax = ax0)

cb.ax.tick_params(labelsize=22)
cb.set_label(label='Relative error [%]',size = 24)

for spine in ax0.spines.values():
    spine.set_linewidth(2)

# Panel b: Correlation coefficient
BrBGcentre  = plt.get_cmap('RdYlBu_r')
cont_rho    = [-1, 0, 0.5,0.6,0.7,0.8,0.9,1]
boundstick  = np.arange(-1,1.2,0.2)
bounds      = [-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]
norm        = mpl.colors.BoundaryNorm(bounds, BrBGcentre.N)

k = 20

ax1 = axs[1]
ax1.set_title(r'(b) Correlation coefficient between time-mean $w_{g}$ and $w_{tot}$ at $\sigma$ '+str(round(isopl[0,k],2))+ ' kg m$^{-3}$',fontsize = 28)
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
grd.xlabel_style = {'size': 22}
grd.ylabel_style = {'size': 22}

ax1.contourf(lon, lat, rho[:,:,k],bounds,transform = ccrs.PlateCarree(),cmap = BrBGcentre, norm = norm,zorder=10)

# Highlight contours
hc = ax1.contour(lon, lat, rho[:,:,k], levels=[0.7, 0.9],colors='black', linewidths=1.5, linestyles=['--', '-'],transform=ccrs.PlateCarree(), zorder=50)

land = cfeature.NaturalEarthFeature('physical', 'land', scale='50m', edgecolor='none', facecolor=cfeature.COLORS['land'], linewidth=.25,zorder = 100)
ax1.add_feature(land, facecolor='k')

cbar = plt.colorbar(
        mpl.cm.ScalarMappable(cmap=BrBGcentre,norm=norm),
    extend='both',
    ticks=boundstick,
    spacing='proportional',
    orientation='vertical',
    shrink=0.8,
    ax = ax1
)

cbar.ax.tick_params(labelsize=22)
cbar.set_label(label='Correlation coefficient',size = 24)

for spine in ax1.spines.values():
    spine.set_linewidth(2)

# Panel c: Vertical gradient proxi

file_b        = sio.loadmat('C:/Users/yago_/Documents/LOCEAN/Codes/vert_grad_w_wglvb_occitens_glob_025_annual_isolevm_filtr.mat')
lon           = file_b['LONimone']
lat           = file_b['LATimone']
wtot_int      = file_b['m_55_tot_int']
wlvb_int      = file_b['m_55_glvb_int']

mask_pos_grd  = np.full(wtot_int.shape, np.nan)
mask_pos_grd[wtot_int>0] = 1
mask_pos_grd[abs(lat)<5] = np.nan

# Bounds:
cmap_LONM          = cm.get_cmap('BrBG_r')
boundticks_rel     = np.arange(-1,1.1,0.1)*100
norm_rel           = mpl.colors.BoundaryNorm(boundticks_rel, cmap_LONM.N)

ee_int             = (wtot_int-wlvb_int)/wtot_int

ax0 = axs[2]
ax0.set_title(r'(c) Relative difference in vertical gradient ($\partial_z w_{tot}$ - $\partial_z w_{g}$)/$\partial_z w_{tot}$',fontsize = 28)
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
grd.xlabel_style = {'size': 22}
grd.ylabel_style = {'size': 22}  

cs = ax0.contourf(lon, lat, ee_int*100,transform=ccrs.PlateCarree(),cmap = cmap_LONM,levels=boundticks_rel,vmin = -100, vmax = 100,zorder=1000,extend='both')

land = cfeature.NaturalEarthFeature('physical', 'land', scale='50m', edgecolor='none', facecolor=cfeature.COLORS['land'], linewidth=.25,zorder = 100)
ax0.add_feature(land, facecolor='k')

# Postive gradient mask
plt.rcParams['hatch.color'] = (0,0,0)
plt.rcParams['hatch.linewidth'] = 1.2
ax0.contourf(lon, lat, mask_pos_grd,[0,1],transform = ccrs.PlateCarree(),hatches=['...'],colors='none',extend='both',zorder=1000)

cb = fig.colorbar(mpl.cm.ScalarMappable(cmap=cmap_LONM, norm=norm_rel),extend='both',ticks=np.arange(-1,1.2,0.2)*100,spacing='proportional',orientation='vertical',label='Relative difference (%)',shrink=0.8,ax = ax0)

cb.ax.tick_params(labelsize=22)
cb.set_label(label='Relative error (%)',size = 24)

for spine in ax0.spines.values():
    spine.set_linewidth(2)

plt.savefig('figure3_VF.png', bbox_inches='tight', dpi=300)

# Figure 9: OGCM comparison with Ekman pumping

proj         = ccrs.Robinson(central_longitude=-60)
BrBGcentre   = plt.get_cmap('RdYlBu_r')
boundstick   = np.arange(-1, 1.2, 0.2)
bounds       = np.arange(-1, 1.05, 0.1)
norm         = mpl.colors.BoundaryNorm(bounds, BrBGcentre.N)
titles = [
    r'(a) R($w_{Ek}$, $w_{tot}$) at $\sigma$ = 26 kg m$^{{-3}}$',
    r'(b) R($w_{{tot}}$, $w_{{Ek}}$) − R($w_{{tot}}$, $w_g$) at $\sigma$ = 26 kg m$^{{-3}}$',
    r'(c) R($w_{{tot}}$, $w_g$) at 100m',
    r'(d) R($w_{{Ek}}$, $w_{{tot}}$) at 100m'
]

# (a) R (Ekman vs w) sigma26
file = sio.loadmat('.../OCCITENS/rho_wek_w_occitens_glob_025_annual_isolevm_filtr.mat')
lon, lat, rho_a, isopl = file['LONimone'], file['LATimone'], file['rho_int'], file['isopl']
k_a = 20

# (b) Diff R (Ekman vs w) - (wg vs w) sigma 26
file1 = sio.loadmat('.../OCCITENS/rho_wek_w_occitens_glob_025_annual_isolevm_filtr.mat')
file2 = sio.loadmat('.../OCCITENS/rho_w_wglvb_occitens_glob_025_annual_isolevm_filtr.mat')
rho_b1, rho_b2 = file1['rho_int'], file2['rho_int']
lon, lat = file1['LONimone'], file1['LATimone']
k_b = 20
diff_b = np.clip(rho_b1[:, :,k_b] - rho_b2[:, :,k_b], -1, 1)

# (c) R (wg vs w) 100 m
file = sio.loadmat('.../OCCITENS/rho_wglvb_w_occitens_glob_025_annual_horlev_filtr.mat')
rho_c = file['rho_int'][:, :]

# (d) R (Ekman vs w) 100 m
file = sio.loadmat('.../OCCITENS/rho_wek_w_occitens_glob_025_annual_horlev_filtr.mat')
rho_d = file['rho_int'][:, :]

# Figure:
fig, axs             = plt.subplots(2, 2, figsize=(22, 12), subplot_kw={"projection": proj})
axs                  = axs.ravel()

datasets             = [rho_a[:, :, k_a],diff_b,rho_c,rho_d]
highlight_levels_all = [[0.7, 0.9],[0],[0.7, 0.9],[0.7, 0.9]]
linestyles_all       = [['--', '-'],['-'],['--', '-'],['--', '-']]

for i, ax in enumerate(axs):
    pcm = ax.contourf(lon, lat, datasets[i], levels=bounds, cmap=BrBGcentre, norm=norm,
                      transform=ccrs.PlateCarree(), zorder=10)
    
    hc = ax.contour(lon, lat, datasets[i], levels=highlight_levels_all[i],
                    colors='black', linewidths=1.5, linestyles=linestyles_all[i],
                    transform=ccrs.PlateCarree(), zorder=50)
    
    ax.add_feature(cfeature.NaturalEarthFeature(
        'physical', 'land', '50m', edgecolor='none', facecolor='k', linewidth=0.25), zorder=100)
    
    gl = ax.gridlines(draw_labels=True, xlocs=range(-180, 181, 90), ylocs=range(-60, 61, 30), color='gray', zorder=1000)
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'size': 20}
    gl.ylabel_style = {'size': 20}
    
    for spine in ax.spines.values():
        spine.set_linewidth(2)
    
    ax.set_title(titles[i], fontsize = 24)

# Colorbar
cbar_ax = fig.add_axes([0.92, 0.2, 0.015, 0.6])
cbar = fig.colorbar(mpl.cm.ScalarMappable(cmap=BrBGcentre, norm=norm), cax=cbar_ax,
                    orientation='vertical', ticks=boundstick, extend='both')
cbar.ax.tick_params(labelsize=20)
cbar.set_label('Correlation coefficient', size=22)

plt.subplots_adjust(wspace=0.1, hspace=0.2)

plt.savefig('figure9_VF.png', bbox_inches='tight', dpi=300)


