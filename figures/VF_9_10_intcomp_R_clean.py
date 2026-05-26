import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.io as sio
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.colors import ListedColormap

# Load variables:
meanm   = sio.loadmat('C:\\Users\\yago_\\Documents\\LOCEAN\\OLIV3 paper\\V4 codes\\intercomparison_metrics\\mean_w.mat')

# Grid
lonbox5 = meanm['lon_box5m']
latbox5 = meanm['lat_box5m']
    
# rho and matsig
rho     = sio.loadmat('C:\\Users\\yago_\\Documents\\LOCEAN\\OLIV3 paper\\V4 codes\\intercomparison_metrics\\R_sign_RMSE_w.mat')['rho']
sig     = sio.loadmat('C:\\Users\\yago_\\Documents\\LOCEAN\\OLIV3 paper\\V4 codes\\intercomparison_metrics\\R_sign_RMSE_w.mat')['matsigcorr']

# MLD mask 
mldmask_occ = sio.loadmat('C:\\Users\\yago_\\Documents\\LOCEAN\\Data\\OCCITENS\\occitens_glob_5_mldmask.mat')['mld_cont_bm']
mldmask_oli = sio.loadmat('C:\\Users\\yago_\\Documents\\LOCEAN\\Data\\OCCITENS\\oliv3_glob_5_mldmask.mat')['mld_cont_bm']
mldmask_ecc = sio.loadmat('C:\\Users\\yago_\\Documents\\LOCEAN\\Data\\OCCITENS\\ecco_glob_5_mldmask.mat')['mld_cont_bm'] 

# %% Figure 9

fig, axes_r = plt.subplots(nrows=3, ncols=2,sharex=True,figsize=(18, 14),subplot_kw={"projection": ccrs.Miller(central_longitude=-60)}, layout='constrained',)
proj = ccrs.Robinson(central_longitude=-60)

fig.suptitle('$w$ correlation coefficient at $\sigma$ 26 kg m$^{-3}$', fontsize = 28, y = 1.03)

k = 20

label_text  = ['(a)','(b)','(c)','(d)','(e)','(f)']
title_name  = ['(a) $GLORYS$ $vs$ $ECCO$',       '(b) $OLIV3$ $vs$ $OMEGA3D$',
              '(c) $OLIV3$ $vs$ $ECCO$',        '(d) $OMEGA3D$ $vs$ $ECCO$',
              '(e) $OLIV3$ $vs$ $GLORYS$',      '(f) $OMEGA3D$ $vs$ $GLORYS$',
              ]
iii         = [5,1,5,5,3,3]
jjj         = [3,0,0,1,0,1]

corrm   = sio.loadmat('C:/Users/yago_/Documents/LOCEAN/OLIV3 paper/corr_RMSE_GLVBw3D_w.mat')
RYB2sel = corrm['RYB2_sel']

cmp = ListedColormap(RYB2sel)

cont_rho    = [-1, 0, 0.5,0.6,0.7,0.8,0.9,1]
boundstick  = np.arange(-1,1.2,0.2)
bounds      = [-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]
norm        = mpl.colors.BoundaryNorm(bounds, cmp.N)

for ij, ax0 in enumerate(axes_r.flat):
        
    ax0.set_title(title_name[ij],fontsize=26)
    ax0.set_xlabel(r'$Longitude(^o)$',fontsize=24)
    ax0.set_ylabel(r'$Latitude(^o)$',fontsize=24)
    title_text = title_name[ij]
    Rvar = rho[:,:,k,jjj[ij],iii[ij]]
    
    Rvar[latbox5==0] = np.nan
    Rvar[latbox5==-5] = np.nan
    
    ax0.set_extent([-180,180,-60,60],crs=ccrs.PlateCarree())
    
    ax0.pcolor(lonbox5+2.5,latbox5+2.5, Rvar, transform=ccrs.PlateCarree(), cmap = cmp, norm = norm)
    
    sigk = sig[:,:,k,jjj[ij],iii[ij]]
    sigk[latbox5==0] = np.nan
    sigk[latbox5==-5] = np.nan
    ax0.scatter(lonbox5*sigk+2.5,latbox5*sigk+2.5, s=5, c='black', marker='.', transform=ccrs.PlateCarree(),zorder=1000)

    # MLD mask --------------------
    plt.rcParams['hatch.color'] = 'black'
    plt.rcParams['hatch.linewidth'] = 2.0
    r_mask = np.full(sigk.shape, np.nan)
    r_mask[rho[:,:,k,2,4] < 0.5] = 1
    r_mask[latbox5==0] = np.nan
    r_mask[latbox5==-5] = np.nan
    if ij == 1 or ij == 2 or ij == 4:
        hatch_pcolor = ax0.pcolor(lonbox5 + 2.5, latbox5 + 2.5, r_mask, transform=ccrs.PlateCarree(), cmap="Grays", alpha=0, hatch="\\\\")

    # Coastlines ------------------
    land = cfeature.NaturalEarthFeature('physical', 'land', scale='50m', edgecolor='none', facecolor=cfeature.COLORS['land'], linewidth=.25)
    ax0.add_feature(land, facecolor='k',zorder=10000)
    
    # Grid ------------------------
    gl = ax0.gridlines(draw_labels=True, 
                       xlocs=range(-180, 181, 90), 
                       ylocs=range(-60, 61, 30), 
                       color='gray', zorder=1)

    gl.right_labels     = False if ij % 2 == 0 else True  
    gl.left_labels      = True if ij % 2 == 0 else False  
    gl.bottom_labels    = True if ij >= 4 else False  
    gl.top_labels       = False                            
    gl.xlabel_style     = {"size": 18, "color": "black"}
    gl.ylabel_style     = {"size": 18, "color": "black"}
    
    for spine in ax0.spines.values():
        spine.set_linewidth(2)
    ax0.tick_params(axis="both", width=2, length=6, labelsize=12)

# Colorbar:
cbar = plt.colorbar(
        mpl.cm.ScalarMappable(norm = norm,cmap=cmp),
    extend='both',
    ticks=boundstick,
    spacing='proportional',
    shrink=0.5,
    ax=axes_r,
    location='bottom',
)

cbar.ax.tick_params(labelsize=20)
cbar.set_label('Correlation coefficient', fontsize=22)
#plt.savefig('C:\\Users\\yago_\\Documents\\LOCEAN\\OLIV3 paper\\V4 codes\\figuresVF\\figure7_VF.png', bbox_inches='tight', dpi=300)

# %% Figure 10

depth_levels    = [18, 20, 28]
depth_labels    = [r'$\sigma$ 25.5 kg m$^{-3}$', r'$\sigma$ 26 kg m$^{-3}$', r'$\sigma$ 27 kg m$^{-3}$']
linestyle       = [ ':', '-',  '-','--', '--','--', ]
title_name      = ['$GLORYS$ $vs$ $ECCO$',       
              '$OLIV3$ $vs$ $ECCO$', 
              '$OLIV3$ $vs$ $GLORYS$',  
              '$OLIV3$ $vs$ $OMEGA3D$',
              '$OMEGA3D$ $vs$ $ECCO$',
              '$OMEGA3D$ $vs$ $GLORYS$',
              ]
iii             = [5,5,3,1,5,3]
jjj             = [3,0,0,0,1,1]

# Colormap
num_colors  = 6  
cmap        = plt.cm.turbo 
colors_list = [cmap(i / (num_colors - 1)) for i in range(num_colors)]

legend_handles = []
legend_labels = []

# Plot
fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(7, 9), sharey=True)

for idx_ax, (ax, k, label) in enumerate(zip(axes, depth_levels, depth_labels)):
    for idx, (i, j) in enumerate(zip(iii, jjj)):
        ROCC = rho[:,:,k,2,4]
        mask_R_OCC = np.full(ROCC.shape,np.nan)
        mask_R_OCC[ROCC>0.5] = 1
        rho_k = rho[:, :, k, j, i]
        rho_k[latbox5 == 0] = np.nan
        rho_k[latbox5 == -5] = np.nan
        lat_vals = latbox5[0,:]
        median_rho = np.nanmedian(rho_k*mask_R_OCC, axis=0)

        line, = ax.plot(lat_vals+2.5, median_rho, linestyle = linestyle[idx], label=title_name[idx], color=colors_list[idx], linewidth=2)
        
        if idx_ax == 0:
            legend_handles.append(line)
            legend_labels.append(title_name[idx])

    ax.grid(True)
    
    if idx_ax < 2:  
        ax.set_xticklabels([])
        ax.set_ylabel("")
    else:
        ax.set_xlabel('Latitude', fontsize=14)
    ax.set_xlim([-80,80])
    ax.set_ylim([-0.5,1])
    
    ax.hlines(0,-80,80,color='dimgray',linewidth = 1)
    
    ax.annotate(label,
            xy=(0.03, 0.15),
            xycoords='axes fraction',
            fontsize=13,
            bbox=dict(boxstyle="round,pad=0.4", fc="white",alpha = 0.9, ec="lightgray", lw=1))
    
    for spine in ax.spines.values():
        spine.set_linewidth(2)
    ax.tick_params(axis='both', width=2, length=6, labelsize=12)

axes[1].set_ylabel('Median correlation coefficient', fontsize=14)

fig.legend(handles=legend_handles,
           labels=legend_labels,
           loc='lower center',
           ncol=2,
           fontsize=12,
           title="Dataset Pairs",
           title_fontsize=14,
           frameon=True, 
           fancybox=True, 
           edgecolor='gray',
           bbox_to_anchor=(0.55, -0.15))

plt.suptitle('Median correlation coefficient across various depths', fontsize=18, y=0.98)
plt.tight_layout()