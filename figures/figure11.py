# Correlation coefficient Ekman pumping, w_tot and w_g

import scipy.io as sio
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import cartopy.crs as ccrs
import cartopy.feature as cfeature

proj = ccrs.Robinson(central_longitude=-60)
BrBGcentre = plt.get_cmap('RdYlBu_r')
boundstick = np.arange(-1, 1.2, 0.2)
bounds = np.arange(-1, 1.05, 0.1)
norm = mpl.colors.BoundaryNorm(bounds, BrBGcentre.N)
titles = [
    r'(a) R($w_{Ek}$, $w_{tot}$) at $\sigma$ = 26 kg m$^{{-3}}$',
    r'(b) R($w_{{tot}}$, $w_{{Ek}}$) − R($w_{{tot}}$, $w_g$) at $\sigma$ = 26 kg m$^{{-3}}$',
    r'(c) R($w_{{tot}}$, $w_g$) at 100m',
    r'(d) R($w_{{Ek}}$, $w_{{tot}}$) at 100m'
]

# %% Load data
# (a)
file = sio.loadmat('C:/Users/yago_/Documents/LOCEAN/Codes/rho_wek_w_occitens_glob_025_annual_isolevm_filtr.mat')
lon, lat, rho_a, isopl = file['LONimone'], file['LATimone'], file['rho_int'], file['isopl']
k_a = 20

# (b)
file1 = sio.loadmat('C:/Users/yago_/Documents/LOCEAN/Codes/rho_wek_w_occitens_glob_025_annual_isolevm_filtr.mat')
file2 = sio.loadmat('C:/Users/yago_/Documents/LOCEAN/Codes/rho_w_wglvb_occitens_glob_025_annual_isolevm_filtr.mat')
rho_b1, rho_b2 = file1['rho_int'], file2['rho_int']
lon, lat = file1['LONimone'], file1['LATimone']
k_b = 20
diff_b = np.clip(rho_b1[:, :,k_b] - rho_b2[:, :,k_b], -1, 1)

# (c)
file = sio.loadmat('C:/Users/yago_/Documents/LOCEAN/Codes/rho_wglvb_w_occitens_glob_025_annual_horlev_filtr.mat')
rho_c = file['rho_int'][:, :]

# (d)
file = sio.loadmat('C:/Users/yago_/Documents/LOCEAN/Codes/rho_wek_w_occitens_glob_025_annual_horlev_filtr.mat')
rho_d = file['rho_int'][:, :]

# %% Figure

fig, axs = plt.subplots(2, 2, figsize=(22, 12), subplot_kw={"projection": proj})
axs = axs.ravel()

datasets              = [rho_a[:, :, k_a], diff_b, rho_c, rho_d]
highlight_levels_all  = [[0.7, 0.9], [0], [0.7, 0.9], [0.7, 0.9]]
linestyles_all        = [['--', '-'], ['-'], ['--', '-'], ['--', '-']]

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

plt.savefig('figure9_VF.png',bbox_inches='tight', dpi=300)
