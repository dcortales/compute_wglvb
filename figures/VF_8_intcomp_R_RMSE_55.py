
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
import matplotlib.colors as mcolors
from matplotlib.legend_handler import HandlerLine2D, HandlerTuple
from matplotlib.patches import Rectangle
import matplotlib.patches as mpatches
import matplotlib.colors as colors
from scipy.ndimage.interpolation import zoom
import xycmap
import pandas as pd

# GRID:
meanm   = sio.loadmat('C:\\Users\\yago_\\Documents\\LOCEAN\\OLIV3 paper\\V4 codes\\intercomparison_metrics\\mean_w.mat')
lonbox5 = meanm['lon_box5m']
latbox5 = meanm['lat_box5m']

wm    = meanm['wm']
    
# rho,rmse, and matsig
rho    = sio.loadmat('C:\\Users\\yago_\\Documents\\LOCEAN\\OLIV3 paper\\V4 codes\\intercomparison_metrics\\R_sign_RMSE_w.mat')['rho']
sig    = sio.loadmat('C:\\Users\\yago_\\Documents\\LOCEAN\\OLIV3 paper\\V4 codes\\intercomparison_metrics\\R_sign_RMSE_w.mat')['matsigcorr']
rmse   = sio.loadmat('C:\\Users\\yago_\\Documents\\LOCEAN\\OLIV3 paper\\V4 codes\\intercomparison_metrics\\R_sign_RMSE_w.mat')['RMSE']

# MLD mask 
mldmask_occ = sio.loadmat('C:\\Users\\yago_\\Documents\\LOCEAN\\Data\\OCCITENS\\occitens_glob_5_mldmask.mat')['mld_cont_bm']
mldmask_oli = sio.loadmat('C:\\Users\\yago_\\Documents\\LOCEAN\\Data\\OCCITENS\\oliv3_glob_5_mldmask.mat')['mld_cont_bm']
mldmask_ecc = sio.loadmat('C:\\Users\\yago_\\Documents\\LOCEAN\\Data\\OCCITENS\\ecco_glob_5_mldmask.mat')['mld_cont_bm'] 

lvbmask_occ = sio.loadmat('C:\\Users\\yago_\\Documents\\LOCEAN\\Data\\OCCITENS\\occitens_glob_5_lvbmask.mat')['lvb_bm']


# Bivariate colormap ------------------------------------------------------

corner_colors = ([0.19830834294502114, 0.23114186851211074, 0.5938485198000769, 1.0],
                 [0.6470588235294118, 0.0, 0.14901960784313725, 1.0],
                 [0.008304498269896196, 0.25266051518646673, 0.11449442522106883, 1.0],
                 [0.32941176470588235, 0.18823529411764706, 0.0196078431372549, 1.0]
                 )

n = (10,10)
xn, yn = n
if xn < 2 or yn < 2:
    raise ValueError("Expected n >= 2 categories.")

color_array = np.array(
    [
        [
            list(colors.to_rgba(corner_colors[0])),
            list(colors.to_rgba('ivory')),
            list(colors.to_rgba(corner_colors[1])),
        ],
                 
        [
            list(colors.to_rgba(corner_colors[2])),
            list(colors.to_rgba('gold')),
            list(colors.to_rgba(corner_colors[3])),
        ],
    ],
)
zoom_factor_x = xn / 3  # Divide by the original two categories.
zoom_factor_y = yn / 2
zcolors = zoom(color_array, (zoom_factor_y, zoom_factor_x, 1), order=1)

sx = pd.Series(np.arange(-1,1.1,0.1))
sx = pd.to_numeric(sx)
sy = pd.Series(np.arange(0,10.5,0.5))
sy = pd.to_numeric(sy)


cmap2 = xycmap.custom_xycmap(corner_colors = corner_colors,n=n)
fig = plt.figure(figsize=(7,7))
cax = fig.add_axes([1,0.25,0.5,0.5])
cax = xycmap.bivariate_legend(ax=cax, sx=sx, sy=sy, cmap=cmap2)
zcolors = cmap2
zcolors[zcolors>1] = 1
zcolors[zcolors<0] = 0
cax.set_xlabel(r'$Correlation$ $coeffcient$',size=12)
cax.set_ylabel(r'$RMSE$ $[10^{-7}m$ $s^{-1)}]$',size=12)
#plt.savefig('R_RMSE_colorbar.png', bbox_inches='tight', dpi=300)

# Limits
Rb      = np.arange(-1,1.2,0.2)
RMSEb   = np.arange(0,11,1)

# ----------------------------------------------------------------------

fig, axes_r     = plt.subplots(nrows=4, ncols=3,sharex=True,figsize=(25, 17),subplot_kw={"projection": ccrs.Miller(central_longitude=-60)}, layout='constrained',)
proj            = ccrs.Robinson(central_longitude=-60)

fig.suptitle('$w$ correlation coefficient and RMSE at $\sigma$ 26 kg m$^{-3}$', fontsize = 28, y = 1.05)

k = 20

label_text = ['(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)','(i)','(j)','(i)']
title_name = ['$OLIV3$ $vs$ $NEMO$',        '$OMEGA3D$ $vs$ $NEMO$',        '$OLIV3$ $vs$ $OMEGA3D$',
              '$OLIV3$ $vs$ $NEMO$ $w_g$',  '$OMEGA3D$ $vs$ $NEMO$ $w_g$',  '$GLORYS$ $vs$ $ECCO$',
              '$OLIV3$ $vs$ $ECCO$',        '$OMEGA3D$ $vs$ $ECCO$', 'nan',
              '$OLIV3$ $vs$ $GLORYS$',      '$OMEGA3D$ $vs$ $GLORYS$',
              ]

# iii = [1,4,4,2,4,5,3,4,0,5,5]
# jjj = [0,1,0,0,2,3,0,3,0,0,4]

iii = [2,2,1,4,4,5,5,5,0,3,3]
jjj = [0,1,0,0,1,3,0,1,0,0,1]


bounds_RMSE     = [0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10]
cmap1           = cm.get_cmap('viridis')

norm3 = mpl.colors.BoundaryNorm(bounds_RMSE, cmap1.N)

for ij, ax0 in enumerate(axes_r.flat):
    
    # Field selection: 
        
    if ij == 8:
        ax0.axis('off')
    elif ij == 11:
        ax0.axis('off')
    else:
        ax0.set_title(title_name[ij],fontsize=26)
        ax0.set_xlabel(r'$Longitude(^o)$',fontsize=24)
        ax0.set_ylabel(r'$Latitude(^o)$',fontsize=24)
        title_text = title_name[ij]
        Rvar, RMSEvar = rho[:,:,k,jjj[ij],iii[ij]], rmse[:,:,k,iii[ij],jjj[ij]]*10**(7)
        
        Rvar[latbox5==0] = np.nan
        Rvar[latbox5==-5] = np.nan
        
        RMSEvar[latbox5==0] = np.nan
        RMSEvar[latbox5==-5] = np.nan
        
        ax0.set_extent([-180,180,-60,60],crs=ccrs.PlateCarree())
        
        bcmp_ind = np.zeros((Rvar.shape[0],Rvar.shape[1],2))*np.nan
        for ii in np.arange(0,Rvar.shape[0]):
            for jj in np.arange(0,Rvar.shape[1]):
                if np.isnan(Rvar[ii,jj]) == False:
                    for rr in np.arange(0,len(Rb)-1):
                        if Rvar[ii,jj] >= Rb[rr] and Rvar[ii,jj] < Rb[rr+1]:
                            bcmp_ind[ii,jj,0] = rr
                    for mm in np.arange(0,len(RMSEb)-1):
                        if RMSEvar[ii,jj] >= RMSEb[mm] and RMSEvar[ii,jj] < RMSEb[mm+1]:
                            bcmp_ind[ii,jj,1] = mm  
                    if RMSEvar[ii,jj] >= RMSEb[len(RMSEb)-1]:
                        bcmp_ind[ii,jj,1] = len(RMSEb)-2  
        
        ax0.scatter(lonbox5[0,0]+2.5,latbox5[0,0]+2.5,s = 30,color = 'white',marker = "s", transform=ccrs.PlateCarree())
        for ii in np.arange(0,Rvar.shape[0]):
            for jj in np.arange(0,Rvar.shape[1]):
                if np.isnan(bcmp_ind[ii,jj,0]) == False:
                    #ax0.scatter(lonbox5[ii,jj]+2.5,latbox5[ii,jj]+2.5,s = 30,color = zcolors[int(bcmp_ind[ii,jj,1])-1,int(bcmp_ind[ii,jj,0])-1,:],marker = "s", transform=ccrs.PlateCarree())
                    ax0.scatter(lonbox5[ii,jj],latbox5[ii,jj],s = 1,color = 'white',marker = "s", transform=ccrs.PlateCarree())
                    xy = [lonbox5[ii,jj],latbox5[ii,jj]]
                    ax0.add_patch(mpatches.Rectangle(xy=[lonbox5[ii,jj],latbox5[ii,jj]],
                                                      width=5,height=5,
                                                      facecolor = zcolors[int(bcmp_ind[ii,jj,1]),int(bcmp_ind[ii,jj,0])-1,:],
                                                      edgecolor = 'none', 
                                                      transform=ccrs.PlateCarree(),
                                                      zorder = 1000))
        sigk = sig[:,:,k,jjj[ij],iii[ij]]
        sigk[latbox5==0] = np.nan
        sigk[latbox5==-5] = np.nan
        ax0.scatter(lonbox5*sigk+2.5,latbox5*sigk+2.5, s=5, c='black', marker='.', transform=ccrs.PlateCarree(),zorder=1000)

    # MLD mask --------------------
        plt.rcParams['hatch.color'] = 'white'  # Choose your hatch color here
        plt.rcParams['hatch.linewidth'] = 2.0
        if ij == 0 or ij == 1 or ij == 3 or ij == 4 or ij == 9 or ij == 10:
            combined_mask = np.where(np.isnan(mldmask_occ[:, :, k]) & np.isnan(mldmask_oli[:, :, k]), np.nan, 1)
            hatch_pcolor = ax0.pcolor(
               lonbox5 + 2.5, latbox5 + 2.5, combined_mask, 
               transform=ccrs.PlateCarree(), cmap="Grays", alpha=0, hatch="///",
               zorder = 10000
           )
        if ij == 2:
            hatch_pcolor = ax0.pcolor(
               lonbox5 + 2.5, latbox5 + 2.5, mldmask_oli[:, :, k], 
               transform=ccrs.PlateCarree(), cmap="Grays", alpha=0, hatch="///",
               zorder = 10000
           )
        if ij == 6 or ij == 7:
            combined_mask = np.where(np.isnan(mldmask_ecc[:, :, k]) & np.isnan(mldmask_oli[:, :, k]), np.nan, 1)
            hatch_pcolor = ax0.pcolor(
               lonbox5 + 2.5, latbox5 + 2.5, combined_mask, 
               transform=ccrs.PlateCarree(), cmap="Grays", alpha=0, hatch="///",
               zorder = 10000
           )
        if ij == 5:
            combined_mask = np.where(np.isnan(mldmask_occ[:, :, k]) & np.isnan(mldmask_ecc[:, :, k]), np.nan, 1)
            hatch_pcolor = ax0.pcolor(
               lonbox5 + 2.5, latbox5 + 2.5, combined_mask, 
               transform=ccrs.PlateCarree(), cmap="Grays", alpha=0, hatch="///",
               zorder = 10000
           )

    # Coastlines ------------------
        land = cfeature.NaturalEarthFeature('physical', 'land', scale='50m', edgecolor='none', facecolor=cfeature.COLORS['land'], linewidth=.25)
        ax0.add_feature(land, facecolor='k',zorder=10000)
    
    # Label -----------------------
        props = dict(boxstyle='round', facecolor='white', edgecolor='white', alpha=1)  # bbox features
        ax0.text(0, 1.13, label_text[ij], transform=ax0.transAxes, fontsize=26, verticalalignment='top', bbox=props)
    # Grid ------------------------

        gl = ax0.gridlines(draw_labels=True, 
                   xlocs=range(-180, 181, 90), 
                   ylocs=range(-60, 61, 30), 
                   color='gray', zorder=1)

        # Only first column gets left labels (ii % 3 == 0)
        gl.left_labels = (ij % 3 == 0)
        
        # Only last row gets bottom labels (ii // 3 == 3)
        gl.bottom_labels = (ij // 3 == 3)
        
        # No top or right labels
        gl.top_labels = False
        gl.right_labels = False
        
        # Label style
        gl.xlabel_style = {"size": 16, "color": "black"}
        gl.ylabel_style = {"size": 16, "color": "black"}
        
        for spine in ax0.spines.values():
            spine.set_linewidth(2)
        ax0.tick_params(axis="both", width=2, length=6, labelsize=20)

cmap2 = xycmap.custom_xycmap(corner_colors = corner_colors,n=n)
cax = fig.add_axes([0.7,0.1,0.3,0.3])
cax = xycmap.bivariate_legend(ax=cax, sx=sx, sy=sy, cmap=cmap2)
zcolors = cmap2
zcolors[zcolors>1] = 1
zcolors[zcolors<0] = 0
cax.set_xlabel(r'$Correlation$ $coeffcient$',size=20)
cax.set_ylabel(r'$RMSE$ $[10^{-7}m$ $s^{-1)}]$',size=20)
cax.set_yticklabels(np.arange(0,11))
for spine in cax.spines.values():
    spine.set_linewidth(2)
cax.tick_params(axis="both", width=2, length=6, labelsize=16)

plt.savefig('8_R_RMSE_intercomp_55_V4.png', bbox_inches='tight', dpi=300)



## pdf from histogram -----------------------------------------------------
import copy
import seaborn as sns

# Reshape
bins_lim    = [-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]
k           = 25
width       = 0.25
R_OME_OLIr = np.reshape(rho[:,:,k,0,1],(73*33,1))
R_OCC_OLIr = np.reshape(rho[:,:,k,0,2],(73*33,1))
R_GLO_OLIr = np.reshape(rho[:,:,k,0,3],(73*33,1))
R_IDI_OLIr = np.reshape(rho[:,:,k,0,4],(73*33,1))
R_ECC_OLIr = np.reshape(rho[:,:,k,0,5],(73*33,1))

R_GLO_ECCr = np.reshape(rho[:,:,k,3,5],(73*33,1))

R_OCC_OMEr = np.reshape(rho[:,:,k,1,2],(73*33,1))
R_GLO_OMEr = np.reshape(rho[:,:,k,1,3],(73*33,1))
R_IDI_OMEr = np.reshape(rho[:,:,k,1,4],(73*33,1))
R_ECC_OMEr = np.reshape(rho[:,:,k,1,5],(73*33,1))


viridis = plt.get_cmap('viridis')



ax1 = []
xs,ys = np.zeros([10,200]), np.zeros([10,200])
var = [R_OCC_OLIr, R_IDI_OLIr, R_GLO_OLIr, R_ECC_OLIr, R_OME_OLIr, R_GLO_ECCr, R_OCC_OMEr, R_IDI_OMEr, R_ECC_OMEr, R_GLO_OMEr]
for vv in var:
    vv[latbox5==0] = np.nan
    vv[latbox5==-5] = np.nan

for ii in np.arange(0,len(var)):
    ax = plt.hist(var[ii],bins=np.arange(-1,1.1,0.1))
    ax1.append(ax[0])
    plt.figure()
    kdeline     = copy.deepcopy(sns.kdeplot(var[ii],bw_adjust=0.8).lines[0])
    ys[ii,:]    = copy.deepcopy(kdeline.get_ydata())
    xs[ii,:]    = copy.deepcopy(kdeline.get_xdata())
    
    
pb_ax1 = copy.deepcopy(ax1)
for ii in np.arange(0,len(ax1)):
    axj = copy.deepcopy(ax1[ii]/sum(ax1[ii]))
    for jj in np.arange(0,len(axj)):
        if jj == 0:
            pb_ax1[ii][jj] = axj[jj]
        elif jj == len(ax1[ii])-1: 
            pb_ax1[ii][jj] = axj[jj]
        else:
            pb_ax1[ii][jj] = sum(axj[jj-1:jj+2])/3    
    
newcolors = viridis(np.linspace(0,0.99,len(ax1)))

plt.figure(figsize=(15,15))
for ii in np.arange(0,len(xs)):
    plt.plot(np.arange(-1,1,0.1)+0.05,ax1[ii]/sum(ax1[ii]),color=newcolors[ii,:])
    plt.plot(xs[ii,:],ys[ii,:]/10,'--',color=newcolors[ii,:])
    plt.plot(np.arange(-1,1,0.1)+0.05,pb_ax1[ii],color=newcolors[ii,:])
    

viridis = plt.get_cmap('viridis')
newcolors = viridis(np.linspace(0,0.99,8))
            

line_style = ['-', (0,(3,1)), (0,(2,1,1,1)), ':', '-', '-', '-', (0,(3,1)), ':', (0,(2,1,1,1))]
line_color = [newcolors[2,:],newcolors[2,:],newcolors[2,:],newcolors[2,:],newcolors[7,:],'black',newcolors[5,:],newcolors[5,:],newcolors[5,:],newcolors[5,:]]
plt.figure(figsize=(7,15))
ax =plt.axes()
mpl.rcParams['legend.frameon'] = 'False'

plt.hlines(0,0,0.16,color='lightgray',zorder=1)

for ii in np.arange(0,len(xs)):
    print(f'q{ii}')
    globals()[f'q{ii}'], = plt.plot(ys[ii,:]/10,xs[ii,:],linestyle = line_style[ii],color=line_color[ii],linewidth=4,zorder=100)
    #ax.fill_betweenx(xs[ii,:], 0, ys[ii,:]/10, color=line_color[ii],alpha=0.1,zorder = 99)
plt.ylim([-1,1])
plt.xlim([0,0.16])
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.ylabel('Correlation Coeffcient',fontsize=16)
plt.xlabel('Density',fontsize=16)
plt.title('Probability distribution of correlation coefficients \nin intercomparison test',fontsize=18)
var = [R_OCC_OLIr, R_IDI_OLIr, R_GLO_OLIr, R_ECC_OLIr, R_OME_OLIr, R_GLO_ECCr, R_OCC_OMEr, R_IDI_OMEr, R_ECC_OMEr, R_GLO_OMEr]

legend = plt.legend([(q4),(q3),(q0),(q1),(q2),(q8),(q6),(q7),(q9),(q5)],
                     ['OLIV3 vs OMEGA3D','OLIV3 vs ECCO','OLIV3 vs NEMO',r'OLIV3 vs NEMO w$_g$','OLIV3 vs GLORYS',
                      'OMEGA3D vs ECCO','OMEGA3D vs NEMO',r'OMEGA3D vs NEMO w$_g$','OMEGA3D vs GLORYS',
                      'GLORYS vs ECCO'], prop=dict(size='15'), loc='lower right', numpoints=1,
                         handler_map={tuple: HandlerTuple(ndivide=None)}, handlelength = 2, ncol = 1,bbox_to_anchor=(1.4-0.4,0.01))

#plt.savefig('R_pdf.png', bbox_inches='tight', dpi=300)


bins_lim = [-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]
xs,ys = np.zeros([10,11,200]), np.zeros([10,11,200])
for k in np.arange(18,29):

    R_OME_OLIr = np.reshape(rho[:,:,k,0,1],(73*33,1))
    R_OCC_OLIr = np.reshape(rho[:,:,k,0,2],(73*33,1))
    R_GLO_OLIr = np.reshape(rho[:,:,k,0,3],(73*33,1))
    R_IDI_OLIr = np.reshape(rho[:,:,k,0,4],(73*33,1))
    R_ECC_OLIr = np.reshape(rho[:,:,k,0,5],(73*33,1))

    R_GLO_ECCr = np.reshape(rho[:,:,k,3,5],(73*33,1))

    R_OCC_OMEr = np.reshape(rho[:,:,k,1,2],(73*33,1))
    R_GLO_OMEr = np.reshape(rho[:,:,k,1,3],(73*33,1))
    R_IDI_OMEr = np.reshape(rho[:,:,k,1,4],(73*33,1))
    R_ECC_OMEr = np.reshape(rho[:,:,k,1,5],(73*33,1))
    
    latb = np.reshape(latbox5,(73*33,1))
    
    var = [R_OCC_OLIr, R_IDI_OLIr, R_GLO_OLIr, R_ECC_OLIr, R_OME_OLIr, R_GLO_ECCr, R_OCC_OMEr, R_IDI_OMEr, R_ECC_OMEr, R_GLO_OMEr]
    for vv in var:
        vv[latb==0] = np.nan
        vv[latb==-5] = np.nan
        
    for ii in np.arange(0,len(var)):
        #ax = plt.hist(var[ii],bins=np.arange(-1,1.1,0.1))
        #ax1.append(ax[0])
        plt.figure()
        kdeline     = copy.deepcopy(sns.kdeplot(var[ii],bw_adjust=0.8).lines[0])
        ys[ii,k-18,:]    = copy.deepcopy(kdeline.get_ydata())
        xs[ii,k-18,:]    = copy.deepcopy(kdeline.get_xdata())


line_style = ['-', (0,(3,1)), (0,(2,1,1,1)), ':', '-', '-', '-', (0,(3,1)), ':', (0,(2,1,1,1))]
line_color = [newcolors[2,:],newcolors[2,:],newcolors[2,:],newcolors[2,:],newcolors[7,:],'black',newcolors[5,:],newcolors[5,:],newcolors[5,:],newcolors[5,:]]
plt.figure(figsize=(10,15))
ax =plt.axes()
mpl.rcParams['legend.frameon'] = 'False'

for ii in np.arange(0,xs.shape[0]):
    plt.hlines(0.05*ii,-1,1,color='lightgray',zorder=0)

for ii in np.arange(0,xs.shape[0]):
    for ij in np.arange(0,xs.shape[1]):

        plt.plot(xs[ii,ij,:],ys[ii,ij,:]/10+ii/20,linestyle = '-',color='black',linewidth=2,zorder=1/(ii+1),alpha = 1-ij/11)
    #ax.fill_betweenx(xs[ii,:], 0, ys[ii,:]/10, color=line_color[ii],alpha=0.1,zorder = 99)
for ii in np.arange(0,xs.shape[0]):
        plt.plot(xs[ii,2,:],ys[ii,2,:]/10+ii/20,linestyle = '-',color='royalblue',linewidth=2,zorder=1/(ii+1))

plt.xlim([-1,1])
#plt.xlim([0,0.16])
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.ylabel('Correlation Coeffcient',fontsize=16)
plt.xlabel('Density',fontsize=16)
plt.title('Probability distribution of correlation coefficients \nin intercomparison test',fontsize=18)

import pandas as pd

R_OME_OLIr = np.reshape(rho[:,:,k,0,1],(73*33,1))
R_OCC_OLIr = np.reshape(rho[:,:,k,0,2],(73*33,1))
R_GLO_OLIr = np.reshape(rho[:,:,k,0,3],(73*33,1))
R_IDI_OLIr = np.reshape(rho[:,:,k,0,4],(73*33,1))
R_ECC_OLIr = np.reshape(rho[:,:,k,0,5],(73*33,1))

R_GLO_ECCr = np.reshape(rho[:,:,k,3,5],(73*33,1))

R_OCC_OMEr = np.reshape(rho[:,:,k,1,2],(73*33,1))
R_GLO_OMEr = np.reshape(rho[:,:,k,1,3],(73*33,1))
R_IDI_OMEr = np.reshape(rho[:,:,k,1,4],(73*33,1))
R_ECC_OMEr = np.reshape(rho[:,:,k,1,5],(73*33,1))

label_name = ['OLIV3 \nvs NEMO','OLIV3 \nvs NEMO w$_g$','OLIV3 \nvs GLORYS','OLIV3 \nvs ECCO','OLIV3 vs \nOMEGA3D',
 'GLORYS \nvs ECCO','OMEGA3D \nvs NEMO','OMEGA3D \nvs NEMO w$_g$','OMEGA3D \nvs ECCO','OMEGA3D \nvs GLORYS']

dfi,dfii = [],[]
for ii in np.arange(0,xs.shape[0]):
    df = pd.DataFrame(ys[ii,2,:]/10)
    aa = df.idxmax()
    dfi.append(aa[0])
    dfii.append(xs[ii,2,aa[0]])
print(df.idxmax())


# dfi,dfii = [],[]
# for ii in np.arange(0,xs.shape[0]):
#     df = pd.DataFrame(ys[ii,2,:]/10)
#     #aa = df.idxmax()
#     aa = sum(ys[ii,2,:]/10*xs[ii,2,:])
#     dfi.append(aa)
#     dfii.append(aa)
# print(df.idxmax())


#dfii = copy.deepcopy(xs[ii,2,dfi])
dfind = np.zeros(len(dfi))
for ii in np.arange(0,len(dfi)):
    mx = 0
    for ij in np.arange(0,len(dfi)):
        if dfii[ij] > mx:
            mx = dfii[ij]
            ind = ij
    dfii[ind] = 0
    dfind[ii] = int(ind)

line_style = ['-', (0,(3,1)), (0,(2,1,1,1)), ':', '-', '-', '-', (0,(3,1)), ':', (0,(2,1,1,1))]
line_color = [newcolors[2,:],newcolors[2,:],newcolors[2,:],newcolors[2,:],newcolors[7,:],'black',newcolors[5,:],newcolors[5,:],newcolors[5,:],newcolors[5,:]]
plt.figure(figsize=(10,15))
ax =plt.axes()
mpl.rcParams['legend.frameon'] = 'False'

for ii in np.arange(0,xs.shape[0]):
    plt.hlines(0.05*ii,-1,1,color='lightgray',zorder=0)

jj = -1
for ii in dfind:
    ii = int(ii)
    jj +=1
    for ij in np.arange(0,xs.shape[1]):
        plt.plot(xs[ii,ij,:],ys[ii,ij,:]/10+jj/20,linestyle = '-',color='black',linewidth=2,zorder=1/(jj+1),alpha = 1-ij/11)
    #ax.fill_betweenx(xs[ii,:], 0, ys[ii,:]/10, color=line_color[ii],alpha=0.1,zorder = 99)
jj = -1
for ii in dfind:
    ii = int(ii)
    jj +=1
    plt.plot(xs[ii,2,:],ys[ii,2,:]/10+jj/20,linestyle = '-',color='royalblue',linewidth=2,zorder=1/(jj+1))
    plt.text(-0.97,jj/20 + 0.01,label_name[ii],fontsize=15)

plt.xlim([-1,1])
plt.ylim([0,0.65])
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.xlabel('Correlation Coeffcient',fontsize=16)
plt.ylabel('Density',fontsize=16)
plt.title('Probability distribution of correlation coefficients \nin intercomparison test',fontsize=18)

plt.savefig('R_pdf_dep_55_V4.png', bbox_inches='tight', dpi=300)

# %% Only correlation coefficients


fig, axes_r     = plt.subplots(nrows=4, ncols=3,sharex=True,figsize=(25, 17),subplot_kw={"projection": ccrs.Miller(central_longitude=-60)}, layout='constrained',)
proj            = ccrs.Robinson(central_longitude=-60)

fig.suptitle('$w$ correlation coefficient at $\sigma$ 26 kg m$^{-3}$', fontsize = 28, y = 1.05)

k = 20

label_text = ['(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)','(i)','(j)','(i)']
title_name = ['$OLIV3$ $vs$ $NEMO$',        '$OMEGA3D$ $vs$ $NEMO$',        '$OLIV3$ $vs$ $OMEGA3D$',
              '$OLIV3$ $vs$ $NEMO$ $w_g$',  '$OMEGA3D$ $vs$ $NEMO$ $w_g$',  '$GLORYS$ $vs$ $ECCO$',
              '$OLIV3$ $vs$ $ECCO$',        '$OMEGA3D$ $vs$ $ECCO$', 'nan',
              '$OLIV3$ $vs$ $GLORYS$',      '$OMEGA3D$ $vs$ $GLORYS$',
              ]

# iii = [1,4,4,2,4,5,3,4,0,5,5]
# jjj = [0,1,0,0,2,3,0,3,0,0,4]

iii = [2,2,1,4,4,5,5,5,0,3,3]
jjj = [0,1,0,0,1,3,0,1,0,0,1]


bounds_RMSE     = [0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10]
cmap1           = cm.get_cmap('viridis')

norm3 = mpl.colors.BoundaryNorm(bounds_RMSE, cmap1.N)

corrm   = sio.loadmat('C:/Users/yago_/Documents/LOCEAN/OLIV3 paper/corr_RMSE_GLVBw3D_w.mat')
RYB2sel = corrm['RYB2_sel']

cmp = ListedColormap(RYB2sel)


#cmp         = ListedColormap(cmpcorr)
cont_rho    = [-1, 0, 0.5,0.6,0.7,0.8,0.9,1]
boundstick  = np.arange(-1,1.2,0.2)
bounds      = [-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]
norm        = mpl.colors.BoundaryNorm(bounds, cmp.N)

for ij, ax0 in enumerate(axes_r.flat):
    
    # Field selection: 
        
    if ij == 8:
        ax0.axis('off')
    elif ij == 11:
        ax0.axis('off')
    else:
        ax0.set_title(title_name[ij],fontsize=26)
        ax0.set_xlabel(r'$Longitude(^o)$',fontsize=24)
        ax0.set_ylabel(r'$Latitude(^o)$',fontsize=24)
        title_text = title_name[ij]
        Rvar, RMSEvar = rho[:,:,k,jjj[ij],iii[ij]], rmse[:,:,k,iii[ij],jjj[ij]]*10**(7)
        
        Rvar[latbox5==0] = np.nan
        Rvar[latbox5==-5] = np.nan
        
        RMSEvar[latbox5==0] = np.nan
        RMSEvar[latbox5==-5] = np.nan
        
        ax0.set_extent([-180,180,-60,60],crs=ccrs.PlateCarree())
        
        ax0.pcolor(lonbox5+2.5,latbox5+2.5, Rvar, transform=ccrs.PlateCarree(), cmap = cmp, norm = norm)
        
      
        sigk = sig[:,:,k,jjj[ij],iii[ij]]
        sigk[latbox5==0] = np.nan
        sigk[latbox5==-5] = np.nan
        ax0.scatter(lonbox5*sigk+2.5,latbox5*sigk+2.5, s=5, c='black', marker='.', transform=ccrs.PlateCarree(),zorder=1000)

    # MLD mask --------------------
        plt.rcParams['hatch.color'] = 'white'  # Choose your hatch color here
        plt.rcParams['hatch.linewidth'] = 2.0
        r_mask = np.full(sigk.shape, np.nan)
        r_mask[rho[:,:,k,2,4] < 0.5] = 1 # Mask correlation coefficient nemo vs nemo wg
        hatch_pcolor = ax0.pcolor(
           lonbox5 + 2.5, latbox5 + 2.5, r_mask, 
           transform=ccrs.PlateCarree(), cmap="Grays", alpha=0, hatch="///",
           zorder = 10000
       )
        
        # if ij == 0 or ij == 1 or ij == 3 or ij == 4 or ij == 9 or ij == 10:
        #     combined_mask = np.where(np.isnan(mldmask_occ[:, :, k]) & np.isnan(mldmask_oli[:, :, k]), np.nan, 1)
        #     hatch_pcolor = ax0.pcolor(
        #        lonbox5 + 2.5, latbox5 + 2.5, combined_mask, 
        #        transform=ccrs.PlateCarree(), cmap="Grays", alpha=0, hatch="///",
        #        zorder = 10000
        #    )
        # if ij == 2:
        #     hatch_pcolor = ax0.pcolor(
        #        lonbox5 + 2.5, latbox5 + 2.5, mldmask_oli[:, :, k], 
        #        transform=ccrs.PlateCarree(), cmap="Grays", alpha=0, hatch="///",
        #        zorder = 10000
        #    )
        # if ij == 6 or ij == 7:
        #     combined_mask = np.where(np.isnan(mldmask_ecc[:, :, k]) & np.isnan(mldmask_oli[:, :, k]), np.nan, 1)
        #     hatch_pcolor = ax0.pcolor(
        #        lonbox5 + 2.5, latbox5 + 2.5, combined_mask, 
        #        transform=ccrs.PlateCarree(), cmap="Grays", alpha=0, hatch="///",
        #        zorder = 10000
        #    )
        # if ij == 5:
        #     combined_mask = np.where(np.isnan(mldmask_occ[:, :, k]) & np.isnan(mldmask_ecc[:, :, k]), np.nan, 1)
        #     hatch_pcolor = ax0.pcolor(
        #        lonbox5 + 2.5, latbox5 + 2.5, combined_mask, 
        #        transform=ccrs.PlateCarree(), cmap="Grays", alpha=0, hatch="///",
        #        zorder = 10000
        #    )

    # Coastlines ------------------
        land = cfeature.NaturalEarthFeature('physical', 'land', scale='50m', edgecolor='none', facecolor=cfeature.COLORS['land'], linewidth=.25)
        ax0.add_feature(land, facecolor='k',zorder=10000)
    
    # Label -----------------------
        props = dict(boxstyle='round', facecolor='white', edgecolor='white', alpha=1)  # bbox features
        ax0.text(0, 1.13, label_text[ij], transform=ax0.transAxes, fontsize=26, verticalalignment='top', bbox=props)
    # Grid ------------------------

        gl = ax0.gridlines(draw_labels=True, 
                   xlocs=range(-180, 181, 90), 
                   ylocs=range(-60, 61, 30), 
                   color='gray', zorder=1)

        # Only first column gets left labels (ii % 3 == 0)
        gl.left_labels = (ij % 3 == 0)
        
        # Only last row gets bottom labels (ii // 3 == 3)
        gl.bottom_labels = ((ij // 3 == 3) | (ij == 5))
    
        # No top or right labels
        gl.top_labels = False
        gl.right_labels = False
        
        # Label style
        gl.xlabel_style = {"size": 16, "color": "black"}
        gl.ylabel_style = {"size": 16, "color": "black"}
        
        for spine in ax0.spines.values():
            spine.set_linewidth(2)
        ax0.tick_params(axis="both", width=2, length=6, labelsize=20)


cbar_ax = fig.add_axes([0.7, 0.05, 0.015, 0.4])  # [left, bottom, width, height]
cb = fig.colorbar(mpl.cm.ScalarMappable(norm = norm,cmap=cmp),
                  cax=cbar_ax,
                  orientation='vertical',
                  extend='both',
                  ticks=boundstick,)

cb.ax.tick_params(labelsize=17)
cb.set_label('Correlation coefficient', fontsize=20)
#plt.savefig('8_R_intercomp_55_V4.png', bbox_inches='tight', dpi=300)

# %% No NEMO

fig, axes_r = plt.subplots(nrows=3, ncols=2,sharex=True,figsize=(18, 14),subplot_kw={"projection": ccrs.Miller(central_longitude=-60)}, layout='constrained',)
proj = ccrs.Robinson(central_longitude=-60)

fig.suptitle('$w$ correlation coefficient at $\sigma$ 26 kg m$^{-3}$', fontsize = 28, y = 1.03)

k = 20

label_text = ['(a)','(b)','(c)','(d)','(e)','(f)']
title_name = ['$OLIV3$ $vs$ $NEMO$',        '$OMEGA3D$ $vs$ $NEMO$',        '$OLIV3$ $vs$ $OMEGA3D$',
              '$OLIV3$ $vs$ $NEMO$ $w_g$',  '$OMEGA3D$ $vs$ $NEMO$ $w_g$',  '$GLORYS$ $vs$ $ECCO$',
              '$OLIV3$ $vs$ $ECCO$',        '$OMEGA3D$ $vs$ $ECCO$', 'nan',
              '$OLIV3$ $vs$ $GLORYS$',      '$OMEGA3D$ $vs$ $GLORYS$',
              ]

title_name = ['(a) $GLORYS$ $vs$ $ECCO$',       '(b) $OLIV3$ $vs$ $OMEGA3D$',
              '(c) $OLIV3$ $vs$ $ECCO$',        '(d) $OMEGA3D$ $vs$ $ECCO$',
              '(e) $OLIV3$ $vs$ $GLORYS$',      '(f) $OMEGA3D$ $vs$ $GLORYS$',
              ]

# iii = [1,4,4,2,4,5,3,4,0,5,5]
# jjj = [0,1,0,0,2,3,0,3,0,0,4]

iii = [2,2,1,4,4,5,5,5,0,3,3]
jjj = [0,1,0,0,1,3,0,1,0,0,1]

iii = [5,1,5,5,3,3]
jjj = [3,0,0,1,0,1]




#cmp         = ListedColormap(cmpcorr)
cont_rho    = [-1, 0, 0.5,0.6,0.7,0.8,0.9,1]
boundstick  = np.arange(-1,1.2,0.2)
bounds      = [-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]
norm        = mpl.colors.BoundaryNorm(bounds, cmp.N)

for ij, ax0 in enumerate(axes_r.flat):
    
    # Field selection: 
        
        ax0.set_title(title_name[ij],fontsize=26)
        ax0.set_xlabel(r'$Longitude(^o)$',fontsize=24)
        ax0.set_ylabel(r'$Latitude(^o)$',fontsize=24)
        title_text = title_name[ij]
        Rvar, RMSEvar = rho[:,:,k,jjj[ij],iii[ij]], rmse[:,:,k,iii[ij],jjj[ij]]*10**(7)
        
        Rvar[latbox5==0] = np.nan
        Rvar[latbox5==-5] = np.nan
        
        RMSEvar[latbox5==0] = np.nan
        RMSEvar[latbox5==-5] = np.nan
        
        ax0.set_extent([-180,180,-60,60],crs=ccrs.PlateCarree())
        
        ax0.pcolor(lonbox5+2.5,latbox5+2.5, Rvar, transform=ccrs.PlateCarree(), cmap = cmp, norm = norm)
        
      
        sigk = sig[:,:,k,jjj[ij],iii[ij]]
        sigk[latbox5==0] = np.nan
        sigk[latbox5==-5] = np.nan
        ax0.scatter(lonbox5*sigk+2.5,latbox5*sigk+2.5, s=5, c='black', marker='.', transform=ccrs.PlateCarree(),zorder=1000)

    # MLD mask --------------------
        plt.rcParams['hatch.color'] = 'black'  # Choose your hatch color here
        plt.rcParams['hatch.linewidth'] = 2.0
        r_mask = np.full(sigk.shape, np.nan)
        r_mask[rho[:,:,k,2,4] < 0.5] = 1 # Mask correlation coefficient nemo vs nemo wg
        r_mask[latbox5==0] = np.nan
        r_mask[latbox5==-5] = np.nan
        if ij == 1 or ij == 2 or ij == 4:
            hatch_pcolor = ax0.pcolor(
                lonbox5 + 2.5, latbox5 + 2.5, r_mask, 
                transform=ccrs.PlateCarree(), cmap="Grays", alpha=0, hatch="\\\\"
            )

    # Coastlines ------------------
        land = cfeature.NaturalEarthFeature('physical', 'land', scale='50m', edgecolor='none', facecolor=cfeature.COLORS['land'], linewidth=.25)
        ax0.add_feature(land, facecolor='k',zorder=10000)
    
    # Label -----------------------
        # props = dict(boxstyle='round', facecolor='white', edgecolor='white', alpha=1)  # bbox features
        # ax0.text(0, 1.13, label_text[ij], transform=ax0.transAxes, fontsize=26, verticalalignment='top', bbox=props)
    # Grid ------------------------

    # Grid ------------------------
        gl = ax0.gridlines(draw_labels=True, 
                           xlocs=range(-180, 181, 90), 
                           ylocs=range(-60, 61, 30), 
                           color='gray', zorder=1)

        gl.right_labels = False if ij % 2 == 0 else True  # Right labels only for right column
        gl.left_labels = True if ij % 2 == 0 else False   # Left labels only for left column
        gl.bottom_labels = True if ij >= 4 else False     # Bottom labels for last row
        gl.top_labels = False                             # No top labels to avoid clutter 
        
        gl.xlabel_style = {"size": 18, "color": "black"}
        gl.ylabel_style = {"size": 18, "color": "black"}
        
        for spine in ax0.spines.values():
            spine.set_linewidth(2)
        ax0.tick_params(axis="both", width=2, length=6, labelsize=12)

        # Colorbar:

cbar = plt.colorbar(
        mpl.cm.ScalarMappable(norm = norm,cmap=cmp),
    extend='both',
    ticks=boundstick,
    spacing='proportional',
    #orientation='horizontal',
    shrink=0.5,
    ax=axes_r,
    location='bottom',
    #pad = 0.5,
)


cbar.ax.tick_params(labelsize=20)
cbar.set_label('Correlation coefficient', fontsize=22)
plt.savefig('8_R_intercomp_55_V4_nn.png', bbox_inches='tight', dpi=300)

# %%

fig, axes_r = plt.subplots(nrows=6, ncols=6,sharex=True,figsize=(18, 14),subplot_kw={"projection": ccrs.Miller(central_longitude=-60)}, layout='constrained',)
proj = ccrs.Robinson(central_longitude=-60)

jjj = [0,0,0,0,0,0,1,1,1,1,1,1,2,2,2,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4,5,5,5,5,5,5]
iii = [0,1,2,3,4,5,0,1,2,3,4,5,0,1,2,3,4,5,0,1,2,3,4,5,0,1,2,3,4,5,0,1,2,3,4,5]
for ij, ax0 in enumerate(axes_r.flat):
    Rvar, RMSEvar = rho[:,:,k,jjj[ij],iii[ij]], rmse[:,:,k,iii[ij],jjj[ij]]*10**(7)
    
    Rvar[latbox5==0] = np.nan
    Rvar[latbox5==-5] = np.nan
    
    RMSEvar[latbox5==0] = np.nan
    RMSEvar[latbox5==-5] = np.nan
    
    ax0.set_extent([-180,180,-60,60],crs=ccrs.PlateCarree())
    
    ax0.pcolor(lonbox5+2.5,latbox5+2.5, Rvar, transform=ccrs.PlateCarree(), cmap = cmp, norm = norm)

# %%
fig, axes_r = plt.subplots(nrows=1, ncols=1,sharex=True,figsize=(18, 14),subplot_kw={"projection": ccrs.Miller(central_longitude=-60)}, layout='constrained',)
proj = ccrs.Robinson(central_longitude=-60)
axes_r.pcolor(lonbox5+2.5,latbox5+2.5, ROCC*mask_R_OCC, transform=ccrs.PlateCarree(), cmap = cmp, norm = norm)
# %% Median r as a function of latitude
# Depth indices to show (as in original code)
depth_levels = [18, 20, 28]
depth_labels = [r'$\sigma$ 25.5 kg m$^{-3}$', r'$\sigma$ 26 kg m$^{-3}$', r'$\sigma$ 27 kg m$^{-3}$']

# Dataset pair indices as defined in your original plot
title_name = ['$OLIV3$ $vs$ $OGCM$',        '$OMEGA3D$ $vs$ $OGCM$',        '$OLIV3$ $vs$ $OMEGA3D$',
              '$OLIV3$ $vs$ $OGCM$ $w_g$',  '$OMEGA3D$ $vs$ $OGCM$ $w_g$',  '$GLORYS$ $vs$ $ECCO$',
              '$OLIV3$ $vs$ $ECCO$',        '$OMEGA3D$ $vs$ $ECCO$',
              '$OLIV3$ $vs$ $GLORYS$',      '$OMEGA3D$ $vs$ $GLORYS$',
              ]
title_name = [        '$OLIV3$ $vs$ $OMEGA3D$',
               '$GLORYS$ $vs$ $ECCO$',
              '$OLIV3$ $vs$ $ECCO$',        '$OMEGA3D$ $vs$ $ECCO$',
              '$OLIV3$ $vs$ $GLORYS$',      '$OMEGA3D$ $vs$ $GLORYS$',
              ]
linestyle = ['-','--','-','-','--',':','-','--','-','--']
iii = [1,4,4,2,4,5,3,4,5,5]
jjj = [0,1,0,0,2,3,0,3,0,4]

title_name = [
    '$OLIV3$ vs $OGCM$',
    '$OLIV3$ vs $OGCM$ $w_g$',
    '$OLIV3$ vs $ECCO$',
    '$OLIV3$ vs $GLORYS$',
    '$OLIV3$ vs $OMEGA3D$',  # moved here
    '$OMEGA3D$ vs $OGCM$',
    '$OMEGA3D$ vs $OGCM$ $w_g$',
    '$OMEGA3D$ vs $ECCO$',
    '$OMEGA3D$ vs $GLORYS$',
    '$GLORYS$ vs $ECCO$'
]

title_name = [

    '$OLIV3$ vs $ECCO$',
    '$OLIV3$ vs $GLORYS$',
    '$OLIV3$ vs $OMEGA3D$',  # moved here
    '$OMEGA3D$ vs $ECCO$',
    '$OMEGA3D$ vs $GLORYS$',
    '$GLORYS$ vs $ECCO$'
]

linestyle = ['-', '-', '-', '-', '-', '--', '--', '--', '--', ':']
linestyle = [ ':', '-',  '-','--', '--','--', ]

title_name = ['$GLORYS$ $vs$ $ECCO$',       '$OLIV3$ $vs$ $ECCO$', '$OLIV3$ $vs$ $GLORYS$',  '$OLIV3$ $vs$ $OMEGA3D$',
                     '$OMEGA3D$ $vs$ $ECCO$',
                  '$OMEGA3D$ $vs$ $GLORYS$',
              ]

iii = [5,5,3,1,5,3]
jjj = [3,0,0,0,1,1]


num_colors = 6  # number of dataset pairs
cmap = plt.cm.turbo  # choose colormap

# Sample colors evenly from the colormap range [0,1]
colors_list = [cmap(i / (num_colors - 1)) for i in range(num_colors)]

legend_handles = []
legend_labels = []

# Create the plot
fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(7, 9), sharey=True)



for idx_ax, (ax, k, label) in enumerate(zip(axes, depth_levels, depth_labels)):
    for idx, (i, j) in enumerate(zip(iii, jjj)):
        # Extract rho and compute median over longitude
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
        
    # ax.set_title(f'Median correlation coefficient at {label}', fontsize=14)

    ax.grid(True)
    
    if idx_ax < 2:  # First and second panels
        ax.set_xticklabels([])
        ax.set_ylabel("")
    else:
        ax.set_xlabel('Latitude', fontsize=14)
    ax.set_xlim([-80,80])
    ax.set_ylim([-0.5,1])
    
    ax.hlines(0,-80,80,color='dimgray',linewidth = 1)
    
    ax.annotate(label,
            xy=(0.03, 0.15),  # Position in axes fraction (x=3%, y=92%)
            xycoords='axes fraction',
            fontsize=13,
            bbox=dict(boxstyle="round,pad=0.4", fc="white",alpha = 0.9, ec="lightgray", lw=1))
    
    for spine in ax.spines.values():
        spine.set_linewidth(2)
    ax.tick_params(axis='both', width=2, length=6, labelsize=12)

axes[1].set_ylabel('Median correlation coefficient', fontsize=14)

# Legend outside plot
#fig.legend(title="Dataset Pairs", loc="center right", bbox_to_anchor=(1.15, 0.5), fontsize=12, title_fontsize=13)

fig.legend(handles=legend_handles,
           labels=legend_labels,
           loc='lower center',
           ncol=2,
           fontsize=12,
           title="Dataset Pairs",
           title_fontsize=14,
           frameon=True,  # <- This puts the legend in a box
           fancybox=True,  # Rounded corners
           edgecolor='gray',
           bbox_to_anchor=(0.55, -0.15))

plt.suptitle('Median correlation coefficient across various depths', fontsize=18, y=0.98)
plt.tight_layout()

plt.savefig('R_pdf_dep_55_new_VFng', bbox_inches='tight', dpi=300)

# %%
import numpy as np
area = np.arange(5,105,10)
da = [5,9,11.25,10.5,9.75,9,8.5,8,11,18]
print(sum(da*area/100))

sa = [6.25,11.5,10.5,10,9.75,9.5,7.5,8,10,17]
print(sum(sa*area/100))

ea = [0,4,8,10.5,10.5,9.5,10,11,13,23.5]
print(sum(ea*area/100))

aa = [0,3.5,8,10.5,10,9.25,9.75,12,13,24]
print(sum(aa*area/100))
