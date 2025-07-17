# Perfect Model Test: OCCITENS Annual Means
# This script computes annual means of Ekman pumping and LVB-derived vertical velocities from the OCCITENS dataset.

import numpy as np
import xarray as xr
from ekman_pumping_comp import ekman_pumping_occ, divergence

from netCDF4 import Dataset
from netCDF4 import Dataset as datas
import h5py

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

# %% 
Vg = [] #,Ug,Vg,MLD,WEK,DIVg   = [],[],[],[],[],[],[],[]                            # Complete data matrices


for aa in anno:
    print("Year: " + str(aa))
    wek = np.zeros((1442,1021,12))#,ug,vg,mld,wek,divg = [],[],[],[],[],[],[],[]                          # Annual matrices
    for mm in np.arange(1,13):
        print("Month: " + str(mm))
        # Velocities:
        pathuo      = r'E:\ORCA025.L75-OCCITENS.003\GRID\%s\ORCA025.L75-OCCITENS.003_y%sm%02d.1m_gridU.nc'%(aa,aa,mm)
        pathvo      = r'E:\ORCA025.L75-OCCITENS.003\GRID\%s\ORCA025.L75-OCCITENS.003_y%sm%02d.1m_gridV.nc'%(aa,aa,mm)
        #pathwo      = r'E:\ORCA025.L75-OCCITENS.003\GRID\%s\ORCA025.L75-OCCITENS.003_y%sm%02d.1m_gridW.nc'%(aa,aa,mm)
        
        fileu = xr.open_dataset(pathuo)
        filev = xr.open_dataset(pathvo)
        
        #um   = fileu.variables['vozocrtx'].values.T
        
        txm,tym     = fileu.variables['sozotaux'].values.T,filev.variables['sometauy'].values.T   # Wind stress
        
        # Ekman pumping:
        wekm        = ekman_pumping_occ(txm,tym,ulon,ulat,deptht,at1,at2,bk)
        
        # Geostrophic velocities:
        #pathvgo     = r'E:\ORCA025.L75-OCCITENS.003\GRID\geo\%s\ORCA025.L75-OCCITENS.003_y%sm%02d.1m_gridVgeo.nc'%(aa,aa,mm)
        #pathvgo     = r'E:\ORCA025.L75-OCCITENS.003\GRID\geo\%s\ORCA025.L75-OCCITENS.003_y%sm%02d.1m_gridVgeo.nc'%(aa,aa,mm)
        
        #filevg = xr.open_dataset(pathvgo)
        
        #vgm     = filevg.variables['vomecrty'].values.T
        
        #ugm1,vgm1   = np.zeros([len(ulat[:,0]),len(vlon[0,:]),len(deptht)]), np.zeros([len(ulat[:,0]),len(vlon[0,:]),len(deptht)]), np.zeros([len(ulat[:,0]),len(vlon[0,:]),len(deptht)])
        #for k in np.arange(0,len(deptht)):
        #    print(k)
        #    ugm1[:,:,k] = np.squeeze(ugm[:,:,k].T)
        #    vgm1[:,:,k] = np.squeeze(vgm[:,:,k].T)
        
        # Divergence:
        #divgm       = divergence(vgm,vlon,vlat,deptht,bf,bi,e3t)
        
        # Mixed layer depth:
        #pathmld     = r'F:\ORCA025.L75-OCCITENS.003\GRID\%s/ORCA025.L75-OCCITENS.003_y%sm%02d.1m_gridflxT.nc'%(aa,aa,mm)
        #filemld     = xr.open_dataset(pathmld)
        #mldm        = filemld.variables['somxl010']
        
        #u[:,:,:,mm] = um
        #v[:,:,:,mm] = vm
        wek[:,:,mm-1] = wekm
        
        #ug.append(ugm)
        #vg.append(vgm)
        #divg.append(divgm)
        
        #wek.append(wekm)
        #vg[:,:,:,mm-1] = np.squeeze(vgm)
        
        #mld.append(mldm)
    
    # Save .nc files:
    from netCDF4 import Dataset as datas
    print('Save annual variables:')


    #print('Save period-average variables:')
    #path_out = r'C:\Users\yago_\Documents\LOCEAN\Data\OCCITENS\div_monthly_%s.nc'%(aa)
    #ncfile = datas(path_out,mode='w',format='NETCDF4_CLASSIC')
    #print(path_out)

    # Create dimensions: 
    #lon_dim             = ncfile.createDimension('x', 1442)     # latitude axis
    #lat_dim             = ncfile.createDimension('y', 1021)    # longitude axis
    #dep_dim             = ncfile.createDimension('z', 75)     # depth axis
    #dep_dim             = ncfile.createDimension('m', 12)     # depth axis
           
    # Define two variables with the same names as dimensions,
    # a conventional way to define "coordinate variables".
    #lat_file            = ncfile.createVariable('lat', np.float64, ('x','y',))
    #lat_file.units      = 'degrees_north'
    #lat_file.long_name  = 'latitude'
    #lon_file            = ncfile.createVariable('lon', np.float64, ('x','y',))
    #lon_file.units      = 'degrees_east'
    #lon_file.long_name  = 'longitude'
    #dep_file            = ncfile.createVariable('depth', np.float64, ('z',))
    #dep_file.units      = 'm'
    #dep_file.long_name  = 'depth'
    #mon_file            = ncfile.createVariable('month', np.float64, ('m',))
    #mon_file.units      = 'month'
    #mon_file.long_name  = 'month'
           
    # Define a 3D variable to hold the data
    #u_file              = ncfile.createVariable('u',np.float64,('x','y','z','m'))
    #u_file.units        = 'm/s' # degrees Kelvin
    #u_file.standard_name= 'zonal velocity' 

    #v_file              = ncfile.createVariable('v',np.float64,('x','y','z','m'))
    #v_file.units        = 'm/s' # degrees Kelvin
    #v_file.standard_name= 'meridional velocity' 

    #w_file              = ncfile.createVariable('w',np.float64,('x','y','z','m'))
    #w_file.units        = 'm/s' # degrees Kelvin
    #w_file.standard_name= 'vertical velocity' 

    #ug_file              = ncfile.createVariable('ug',np.float64,('x','y','z','m'))
    #ug_file.units        = 'm/s' # degrees Kelvin
    #ug_file.standard_name= 'geostrophic zonal velocity' 

    #vg_file              = ncfile.createVariable('vg',np.float64,('x','y','z','m'))
    #vg_file.units        = 'm/s' # degrees Kelvin
    #vg_file.standard_name= 'geostrophic meridional velocity' 

    #dg_file              = ncfile.createVariable('div',np.float64,('x','y','z','m'))
    #dg_file.units        = 'm/s' # degrees Kelvin
    #dg_file.standard_name= 'geostrophic velocity divergence integral' 

    #ek_file              = ncfile.createVariable('wek',np.float64,('x','y','m'))
    #ek_file.units        = 'm/s' # degrees Kelvin
    #ek_file.standard_name= 'Ekman pumping vertical velocity' 
           
    # Write latitudes, longitudes.
    # Note: the ":" is necessary in these "write" statements
    #lat_file[:,:],  lon_file[:,:],  dep_file[:],mon_file[:]        = ulat, ulon, deptht,month_num
    #u_file[:,:,:],  v_file[:,:,:],  w_file[:,:,:],      = u,v,w
    #w_file[:,:,:,:]      = w
    #ug_file[:,:,:], vg_file[:,:,:], dg_file[:,:,:]      = Ug, Vg, DIVg
    #dg_file[:,:,:,:]      = divg

    #ncfile.close()   
    
    
    
    #U.append(np.mean(u,2))
    #V.append(np.mean(v,2))
    WEK = np.mean(wek,2)
    
    #Ug.append(ug)
    #Vg.append(vg)
    #DIVg.append(divg)
    
    #WEK.append(wek)

    #MLD.append(mld)
    
     
    # Annual averages:
    # Entire period averages:
    
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
    #u_file              = ncfile.createVariable('u',np.float64,('x','y','z'))
    #u_file.units        = 'm\s' # degrees Kelvin
    #u_file.standard_name= 'zonal velocity' 
    
    #v_file              = ncfile.createVariable('v',np.float64,('x','y','z'))
    #v_file.units        = 'm/s' # degrees Kelvin
    #v_file.standard_name= 'meridional velocity' 
    
    #w_file              = ncfile.createVariable('w',np.float64,('x','y','z'))
    #w_file.units        = 'm/s' # degrees Kelvin
    #w_file.standard_name= 'vertical velocity' 
    
    # ug_file              = ncfile.createVariable('ug',np.float64,('x','y','z'))
    # ug_file.units        = 'm/s' # degrees Kelvin
    # ug_file.standard_name= 'geostrophic zonal velocity' 
    
    # vg_file              = ncfile.createVariable('vg',np.float64,('x','y','z'))
    # vg_file.units        = 'm/s' # degrees Kelvin
    # vg_file.standard_name= 'geostrophic meridional velocity' 
    
    #dg_file              = ncfile.createVariable('div',np.float64,('x','y','z'))
    #dg_file.units        = 'm/s' # degrees Kelvin
    #dg_file.standard_name= 'geostrophic velocity divergence integral' 
    
    ek_file              = ncfile.createVariable('wek',np.float64,('x','y'))
    ek_file.units        = 'm/s' # degrees Kelvin
    ek_file.standard_name= 'Ekman pumping vertical velocity' 
           
    # Write latitudes, longitudes.
    # Note: the ":" is necessary in these "write" statements
    lat_file[:,:],  lon_file[:,:],  dep_file[:]         = ulat, ulon, deptht
    #u_file[:,:,:],  v_file[:,:,:],  w_file[:,:,:],      = U, V, W
    #w_file[:,:,:]      = W
    #ug_file[:,:,:], vg_file[:,:,:], dg_file[:,:,:]      = Ug, Vg, DIVg
    ek_file[:,:]      = WEK
    
    ncfile.close()