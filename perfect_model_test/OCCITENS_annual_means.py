# Perfect Model Test: OCCITENS Annual Means
# This script computes annual means of Ekman pumping and LVB-derived vertical velocities from the OCCITENS dataset.

import numpy as np
import xarray as xr
from ekman_pumping_comp import ekman_pumping_occ, divergence

from netCDF4 import Dataset
from netCDF4 import Dataset as datas

'''
Compute Ekman pumping and w in the ocean interior:
- Data:
    - OCCITENS
- Computation:
    - Ekman pumping
    - divergence of geostrophic flow
'''

# %% Load data

# GRID:
pathu           = r'E:\ORCA025.L75-OCCITENS.003\GRID\1993\ORCA025.L75-OCCITENS.003_y1993m01.1m_gridU.nc' 
pathv           = r'E:\ORCA025.L75-OCCITENS.003\GRID\1993\ORCA025.L75-OCCITENS.003_y1993m01.1m_gridV.nc' 
fileu, filev    = xr.open_dataset(pathu),xr.open_dataset(pathv)
ulon,ulat       = fileu.variables['nav_lon'].values.T, fileu.variables['nav_lat'].values.T
vlon,vlat       = filev.variables['nav_lon'].values.T, filev.variables['nav_lat'].values.T

patht           = r'E:\ORCA025.L75-OCCITENS.003\GRID\1993\ORCA025.L75-OCCITENS.003_y1993m01.1m_gridT.nc' 
filet           = xr.open_dataset(patht)
deptht          = filet.variables['deptht'].values

pathgrid        = r'E:\ORCA025.L75-OCCITENS.003\GRID\ORCA025.L75-OCCITENS_mesh_mask.nc' 
filegrid        = xr.open_dataset(pathgrid)
e1t,e2t,e3t     = filegrid.variables['e1t'].values.T,filegrid.variables['e2t'].values.T,filegrid.variables['e3t'].values.T
e1f,e2f         = filegrid.variables['e1f'].values.T,filegrid.variables['e2f'].values.T
e3w             = filegrid.variables['e3w'].values.T
e1v,e2v,e3v     = filegrid.variables['e1v'].values.T,filegrid.variables['e2v'].values.T,filegrid.variables['e3v'].values.T

filegrid.close()

# Combined scale factors:
bk  = np.zeros((len(e1t[:,0]),len(e1t[0,:]),len(e3t[0,0,:])))
at1 = np.zeros((len(e1t[:,0]),len(e1t[0,:]),len(e3t[0,0,:])))
at2 = np.zeros((len(e1t[:,0]),len(e1t[0,:]),len(e3t[0,0,:])))

bf  = np.zeros((len(e1t[:,0]),len(e1t[0,:]),len(e3t[0,0,:])))
bi  = np.zeros((len(e1t[:,0]),len(e1t[0,:]),len(e3t[0,0,:])))

for k in np.arange(0, len(deptht)):
    bk[:,:,k]   = np.squeeze(e1t)*np.squeeze(e2t)*np.squeeze(e3t[:,:,k])
    at1[:,:,k]  = np.squeeze(e1t)*np.squeeze(e3t[:,:,k])
    at2[:,:,k]  = np.squeeze(e2t)*np.squeeze(e3t[:,:,k])
    bf[:,:,k]   = np.squeeze(e1f)*np.squeeze(e2f)*np.squeeze(e3w[:,:,k])
    bi[:,:,k]   = np.squeeze(e1v)*np.squeeze(e3v[:,:,k])
    
month_num       = np.arange(1,13)

# %% Compute annual means
anno = np.arange(1993,2016)

# %% Geostrophic meridional velocities:

for aa in anno:
    print("Year: " + str(aa))
    vg = np.zeros((1442,1021,75,12)) # Annual matrix
    for mm in np.arange(1,13):
        print("Month: " + str(mm))
        pathvgo     = r'E:\ORCA025.L75-OCCITENS.003\GRID\geo\%s\ORCA025.L75-OCCITENS.003_y%sm%02d.1m_gridVgeo.nc'%(aa,aa,mm)
        filevg      = xr.open_dataset(pathvgo)
        
        vgm         = filevg.variables['vomecrty'].values.T
        vg[:,:,:,mm-1] = vgm
    Vg = np.mean(vg,3)  # Average over months

    # Save .nc file:
    print('Save annual variables:')

    path_out = r'C:\Users\yago_\Documents\LOCEAN\Data\OCCITENS\vgeo_annual_%s.nc'%(aa)
    ncfile = datas(path_out,mode='w',format='NETCDF4_CLASSIC')
    print(path_out)

    # Create dimensions: 
    lon_dim             = ncfile.createDimension('x', 1442)     # latitude axis
    lat_dim             = ncfile.createDimension('y', 1021)    # longitude axis
    dep_dim             = ncfile.createDimension('z', 75)     # depth axis
            
    # Define two variables with the same names as dimensions,
    # a conventional way to define "coordinate variables".
    lat_file            = ncfile.createVariable('lat', np.float64, ('x','y',))
    lat_file.units      = 'degrees_north'
    lat_file.long_name  = 'latitude'
    lon_file            = ncfile.createVariable('lon', np.float64, ('x','y',))
    lon_file.units      = 'degrees_east'
    lon_file.long_name  = 'longitude'
    dep_file            = ncfile.createVariable('depth', np.float64, ('z',))
    dep_file.units      = 'm'
    dep_file.long_name  = 'depth'
            
    # Define a 3D variable to hold the data

    vg_file              = ncfile.createVariable('vg',np.float64,('x','y','z'))
    vg_file.units        = 'm/s' # meters per second
    vg_file.standard_name= 'geostrophic meridional velocity' 
            
    # Write variables
    lat_file[:,:]   = ulat
    lon_file[:,:]   = ulon
    dep_file[:]     = deptht
    vg_file[:,:,:]  = Vg

    ncfile.close()

# %% Total vertical velocity

for aa in anno:
    print("Year: " + str(aa))
    w = np.zeros((1442,1021,75,12)) # Annual matrix
    for mm in np.arange(1,13):
        print("Month: " + str(mm))
        # Velocities:
        pathwo      = r'E:\ORCA025.L75-OCCITENS.003\GRID\%s\ORCA025.L75-OCCITENS.003_y%sm%02d.1m_gridW.nc'%(aa,aa,mm)
        
        filew = xr.open_dataset(pathwo)
        wm         = filevg.variables['vovecrtz'].values.T
        w[:,:,:,mm-1] = wm
    W = np.mean(w,3)  # Average over months

    # Save .nc files:
    print('Save annual variables:')
    
    path_out = r'C:\Users\yago_\Documents\LOCEAN\Data\OCCITENS\w_annual_%s.nc'%(aa)
    ncfile = datas(path_out,mode='w',format='NETCDF4_CLASSIC')
    print(path_out)
    
    # Create dimensions: 
    lon_dim             = ncfile.createDimension('x', 1442)     # latitude axis
    lat_dim             = ncfile.createDimension('y', 1021)    # longitude axis
    dep_dim             = ncfile.createDimension('z', 75)     # depth axis
           
    # Define two variables with the same names as dimensions,
    # a conventional way to define "coordinate variables".
    lat_file            = ncfile.createVariable('lat', np.float64, ('x','y',))
    lat_file.units      = 'degrees_north'
    lat_file.long_name  = 'latitude'
    lon_file            = ncfile.createVariable('lon', np.float64, ('x','y',))
    lon_file.units      = 'degrees_east'
    lon_file.long_name  = 'longitude'
    dep_file            = ncfile.createVariable('depth', np.float64, ('z',))
    dep_file.units      = 'm'
    dep_file.long_name  = 'depth'
           
    # Define a 3D variable to hold the data
    w_file              = ncfile.createVariable('w',np.float64,('x','y','z'))
    w_file.units        = 'm/s' # meters per second
    w_file.standard_name= 'vertical velocity' 
           
    # Write variables:
    lat_file[:,:] = ulat
    lon_file[:,:] = ulon
    dep_file[:]   = deptht  
    w_file[:,:,:] = W
    
    ncfile.close()
# %% Compute Ekman pumping

for aa in anno:
    print("Year: " + str(aa))
    wek = np.zeros((1442,1021,12))                         # Annual matrices
    for mm in np.arange(1,13):
        print("Month: " + str(mm))
        # Velocities:
        pathuo      = r'E:\ORCA025.L75-OCCITENS.003\GRID\%s\ORCA025.L75-OCCITENS.003_y%sm%02d.1m_gridU.nc'%(aa,aa,mm)
        pathvo      = r'E:\ORCA025.L75-OCCITENS.003\GRID\%s\ORCA025.L75-OCCITENS.003_y%sm%02d.1m_gridV.nc'%(aa,aa,mm)
  
        fileu       = xr.open_dataset(pathuo)
        filev       = xr.open_dataset(pathvo)
        
        txm, tym    = fileu.variables['sozotaux'].values.T, filev.variables['sometauy'].values.T   # Wind stress
        
        # Ekman pumping:
        wekm        = ekman_pumping_occ(txm,tym,ulon,ulat,deptht,at1,at2,bk)
        wek[:,:,mm-1] = wekm
    WEK = np.mean(wek,2)
    
    # Save .nc files:
    print('Save annual variables:')

    path_out = r'C:\Users\yago_\Documents\LOCEAN\Data\OCCITENS\wek_annual_%s.nc'%(aa)
    ncfile = datas(path_out,mode='w',format='NETCDF4_CLASSIC')
    print(path_out)

    # Create dimensions: 
    lon_dim             = ncfile.createDimension('x', 1442)     # latitude axis
    lat_dim             = ncfile.createDimension('y', 1021)    # longitude axis
    dep_dim             = ncfile.createDimension('z', 75)     # depth axis
            
    # Define two variables with the same names as dimensions,
    # a conventional way to define "coordinate variables".
    lat_file            = ncfile.createVariable('lat', np.float64, ('x','y',))
    lat_file.units      = 'degrees_north'
    lat_file.long_name  = 'latitude'
    lon_file            = ncfile.createVariable('lon', np.float64, ('x','y',))
    lon_file.units      = 'degrees_east'
    lon_file.long_name  = 'longitude'
    dep_file            = ncfile.createVariable('depth', np.float64, ('z',))
    dep_file.units      = 'm'
    dep_file.long_name  = 'depth'
            
    # Define a 3D variable to hold the data
    ek_file              = ncfile.createVariable('wek',np.float64,('x','y'))
    ek_file.units        = 'm/s' # meters per second
    ek_file.standard_name= 'Ekman pumping vertical velocity' 
            
    # Write variables:
    lat_file[:,:]   = ulat
    lon_file[:,:]   = ulon
    dep_file[:]     = deptht
    ek_file[:,:]    = WEK

    ncfile.close()

# %% Compute geostrophic LVB vertical velocities

for aa in anno:
    print("Year: " + str(aa))
    vg = np.zeros((1442,1021,75))
    # Geostrophic velocities:
    for mm in np.arange(1,13):
        print("Month: " + str(mm))
        
        # Geostrophic velocities:
        pathvgo     = r'E:\ORCA025.L75-OCCITENS.003\GRID\geo\%s\ORCA025.L75-OCCITENS.003_y%sm%02d.1m_gridVgeo.nc'%(aa,aa,mm)
        filevg = xr.open_dataset(pathvgo)
        vgm     = np.squeeze(filevg.variables['vomecrty'].values.T)    
        vg += vgm
    Vg = vg/12
    
    # Ekman pumping:
    pathek      = r'C:\Users\yago_\Documents\LOCEAN\Data\OCCITENS\wek_annual_%d.nc'%(aa)
    fileek      = xr.open_dataset(pathek)
    wekm        = fileek.variables['wek'].values
    
    # Divergence:
    wglvbm      = wglvb(Vg,wekm,vlon,vlat,deptht,bf,bi,e3t)
        
    # Save .nc files:
    print('Save annual variables:')
    
    path_out = r'C:\Users\yago_\Documents\LOCEAN\Data\OCCITENS\wglvb_annual_%s.nc'%(aa)
    ncfile = datas(path_out,mode='w',format='NETCDF4_CLASSIC')
    print(path_out)
    
    # Create dimensions: 
    lon_dim             = ncfile.createDimension('x', 1442)     # latitude axis
    lat_dim             = ncfile.createDimension('y', 1021)    # longitude axis
    dep_dim             = ncfile.createDimension('z', 75)     # depth axis
    
    # Create variables:
    lat_file            = ncfile.createVariable('lat', np.float64, ('x','y',))
    lat_file.units      = 'degrees_north'
    lat_file.long_name  = 'latitude'
    lon_file            = ncfile.createVariable('lon', np.float64, ('x','y',))
    lon_file.units      = 'degrees_east'
    lon_file.long_name  = 'longitude'
    dep_file            = ncfile.createVariable('depth', np.float64, ('z',))
    dep_file.units      = 'm'
    dep_file.long_name  = 'depth'
            
    w_file              = ncfile.createVariable('div',np.float64,('x','y','z'))
    w_file.units        = 'm/s' # meters per second
    w_file.standard_name= 'beta-plane geostrophic vertical velocity' 
            
    lat_file[:,:],  lon_file[:,:],  dep_file[:]     = ulat, ulon, deptht
    w_file[:,:,:]                                   = wglvbm
    
    ncfile.close()
