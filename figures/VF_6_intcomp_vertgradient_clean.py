import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.io as sio
from matplotlib import cm
import cartopy.crs as ccrs
import cartopy.feature as cfeature


# GRID:
meanm   = sio.loadmat('C:\\Users\\yago_\\Documents\\LOCEAN\\OLIV3 paper\\V4 codes\\intercomparison_metrics\\mean_w_womld.mat')
lonbox5 = meanm['lon_box5m']
latbox5 = meanm['lat_box5m']
    
# Mean:
slp     = sio.loadmat('C:\\Users\\yago_\\Documents\\LOCEAN\\OLIV3 paper\\V4 codes\\intercomparison_metrics\\slope_w_55new.mat')['m_55']
slp     = sio.loadmat('C:/Users/yago_/Documents/LOCEAN/intercomparison_metrics\slope_w_55.mat')['m_55']

# %% Figure

label_name = [r'(a) $OLIV3$',r'(b) $OMEGA3D$',r'(c) $GLORYS12v1$',r'(d) $ECCOv4r4$',
              r'(e) |$OMEGA3D$| - $|OLIV3|$',r'(f) $|GLORYS12v1|$ - $|OLIV3|$',r'(g) $|ECCOv4r4|$ - $|OLIV3|$',
              r'(h) $|GLORYS12v1|$ - $|OMEGA3D|$',r'(i) $|ECCOv4r4|$ - $|OMEGA3D|$', r'(j) $|GLORYS12v1|$ - $|ECCOv4r4|$']

fig, axes_r = plt.subplots(nrows=5, ncols=2,sharex=True,figsize=(14, 16),subplot_kw={"projection": ccrs.Miller(central_longitude=-60)}, layout='constrained',)
proj = ccrs.Robinson(central_longitude=-60)

cmap_LONM           = cm.get_cmap('BrBG_r')
boundLONM           = np.arange(-2,2.1,0.2)
norm_LONM           = mpl.colors.BoundaryNorm(boundLONM, cmap_LONM.N)

fig.suptitle(r'Vertical gradient of time-mean $w$ and differences among datasets', fontsize = 24, y = 1.03)
ind = [0,1,3,5,1,3,5,3,5,5]
for ii,ax0 in enumerate(axes_r.flat):
    
    ax0.set_title(label_name[ii],fontsize=22)
    ax0.set_extent([-180,180,-60,60],crs=ccrs.PlateCarree())
    
    if ii < 4:
        wk = slp[:,:,ind[ii]]
    elif ii > 3 and ii < 7:
        wk = abs(slp[:,:,ind[ii]]) - abs(slp[:,:,0])
    elif ii > 6 and ii < 9:
        wk = abs(slp[:,:,ind[ii]]) - abs(slp[:,:,1])
    elif ii > 8:
        wk = abs(slp[:,:,ind[ii]]) - abs(slp[:,:,3])
    
    wk[latbox5==0] = np.nan
    wk[latbox5==-5] = np.nan
    ax0.pcolor(lonbox5+2.5,latbox5+2.5, wk*10**4, transform=ccrs.PlateCarree(), cmap = cmap_LONM, norm = norm_LONM)

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
    
# Colorbar: 
    
boundLONM_ticks     = boundLONM[::4]

cbar = plt.colorbar(
        mpl.cm.ScalarMappable(cmap=cmap_LONM, norm = norm_LONM),
    extend='both',
    ticks=boundLONM_ticks,
    spacing='proportional',
    shrink=0.5,
    ax=axes_r,
    location='bottom',
    pad = 0.03,
)

cbar.ax.tick_params(labelsize=18)
cbar.set_label(label=r'$10^{-4}$ $m$ $day^{-1}$ $m^{-1}$',size=20)

#plt.savefig('C:\\Users\\yago_\\Documents\\LOCEAN\\OLIV3 paper\\V4 codes\\review\\6_vertical_structure_55_VF_diff_corr.png', bbox_inches='tight', dpi=300)
