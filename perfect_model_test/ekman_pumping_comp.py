# Ekman pumping computation
import numpy as np

def ekman_pumping_occ(taux,tauy,lon,lat,deptht,at1,at2,bk):
    
    taux, tauy      = np.squeeze(taux), np.squeeze(tauy)
    
    # Calculamos la w_Ek ----------------------------------------------
    
    # Parameters:
    fcO         = 1.454441043328608e-04         # 2*rotation_rate
    f           = fcO*(np.sin(lat*np.pi/180))   # Coriolis parameter
    rho0        = 1025                          # reference pressure
    f[abs(lat)<3] = np.nan
    
    # Gradiente zonal de tau_y
    #for k in np.arange(0,len(deptht)):
    dtauy       = np.squeeze(at2[1:,:,0])*tauy[1:,:]/f[1:,:] - np.squeeze(at2[:-1,:,0])*tauy[:-1,:]/f[:-1,:]
    #tauyi       = 0.5*(dtauy[:-1,:] + dtauy[1:,:])
    dtauydx     = dtauy/np.squeeze(bk[:-1,:,0])

    # Gradiente meridional de tau_x
    dtaux       = np.squeeze(at1[:,:-1,0])*(taux[:,1:]/f[:,1:]) - np.squeeze(at1[:,:-1,0])*(taux[:,:-1]/f[:,:-1])
    #tauxj       = 0.5*(dtaux[:,:-1]+dtaux[:,1:])
    dtauxdy     = dtaux/np.squeeze(bk[:,:-1,0])

    wek         = (1/rho0)*(np.squeeze(dtauydx[:,:-1])-np.squeeze(dtauxdy[:-1,:]))
    wekm        = np.zeros((1442,1021))
    wekm[:-1,:-1] = wek
    return wekm

def ekman_pumping(taux,tauy,lon,lat):
        
    # Calculamos la w_Ek ----------------------------------------------
    # Parameters:
    fcO         = 1.454441043328608e-04         # 2*rotation_rate
    f           = fcO*(np.sin(lat*np.pi/180))   # Coriolis parameter
    rho0        = 1025                          # reference pressure
    Rt          = 6371229                       # Earth radius
    f[abs(lat)<3] = np.nan
    
    # spatial derivatives:
    dx          = Rt*abs(lon[1:,:]-lon[:-1,:])*np.pi/180;
    dy          = Rt*abs(lat[:,1:]-lat[:,:-1])*np.pi/180;
    
    # Gradiente zonal de tau_y
    dtauy       = tauy[1:,:]/f[1:,:] - tauy[:-1,:]/f[:-1,:]
    dtauydx     = dtauy/dx;
    
    # Gradiente meridional de tau_x
    dtaux       = taux[:,1:]/f[:,1:] - taux[:,:-1]/f[:,:-1]
    dtauxdy     = dtaux/dy;

    wek         = (1/rho0)*(dtauydx[:,:-1]-dtauxdy[:-1,:])
    return wek

def divergence(vg,lon,lat,deptht,bf,bi,e3t):
    
    vg          = np.squeeze(vg)
    # Parameters
    fcO         = 1.454441043328608e-04                                         # 2*rotation_rate
    f           = fcO*(np.sin(lat*np.pi/180))                                   # Coriolis parameter
    f[abs(lat)<3] = np.nan

    # Meridional velocity:
    v1          = bi*vg
    vi          = 0.5*(v1[1:,:,:] + v1[:-1,:,:])                
    vj          = 0.5*(vi[:,1:,:] + vi[:,:-1,:])  
    
    dfdz = []
    # Gradiente zonal de tau_y
    for k in np.arange(0,len(deptht)):
        df          = f[:,1:]-f[:,:-1]
        dfdz.append((df[:-1,:]/bf[:-1,:-1,k]).T)                                   # f gradient:
            
    dfdz = np.array(dfdz).T
    
    betav1      = vj*dfdz
    bb          = 0.5*(betav1[:,1:,:] + betav1[:,:-1,:])                        # betav:
    
    dif_lvb     = np.zeros((1442,1021,75))
    dbbdz       = (bb[:,:,:-1] + bb[:,:,1:])/2
    #dif_zt      = abs(deptht[:-1]-deptht[1:])  
    for k in np.arange(0,len(deptht)-1):
        #dif_zt2D = np.ones((dbbdz.shape[0],dbbdz.shape[1]))*dif_zt
        dif_lvb[:-1,1:-1,k] = ((1/f[:-1,1:-1])*dbbdz[:,:,k])                   # vertical gradient/f
    #dif_lvb = np.array(dif_lvb).T
        
    div         = np.cumsum(dif_lvb*np.squeeze(e3t),axis=2)                              # cumulative sum

    return div

def wglvb(vg,wek,lon,lat,deptht,bf,bi,e3t):
    
    vg          = np.squeeze(vg)
    wek         = np.squeeze(wek)
    # Parameters
    fcO         = 1.454441043328608e-04                                         # 2*rotation_rate
    f           = fcO*(np.sin(lat*np.pi/180))                                   # Coriolis parameter
    f[abs(lat)<3] = np.nan

    # Meridional velocity:
    v1          = bi*vg
    vi          = 0.5*(v1[1:,:,:] + v1[:-1,:,:])                
    vj          = 0.5*(vi[:,1:,:] + vi[:,:-1,:])  
    
    dfdz = []
    # Gradiente zonal de tau_y
    for k in np.arange(0,len(deptht)):
        df          = f[:,1:]-f[:,:-1]
        dfdz.append((df[:-1,:]/bf[:-1,:-1,k]).T)                                   # f gradient:
            
    dfdz = np.array(dfdz).T
    
    betav1      = vj*dfdz
    bb          = 0.5*(betav1[:,1:,:] + betav1[:,:-1,:])                        # betav:
    
    dif_lvb     = np.zeros((1442,1021,74))
    dbbdz       = (bb[:,:,:-1] + bb[:,:,1:])/2
    #dif_zt      = abs(deptht[:-1]-deptht[1:])  
    for k in np.arange(0,len(deptht)-1):
        #dif_zt2D = np.ones((dbbdz.shape[0],dbbdz.shape[1]))*dif_zt
        dif_lvb[:-1,1:-1,k] = -((1/f[:-1,1:-1])*dbbdz[:,:,k])                   # vertical gradient/f
    #dif_lvb = np.array(dif_lvb).T
    dlvb            = np.zeros(vg.shape)
    dlvb[:,:,0]     = wek
    dlvb[:,:,1:]    = dif_lvb*np.squeeze((e3t[:,:,:-1]))
    wglvb           = np.cumsum(dlvb,axis=2)                         # cumulative sum

    return wglvb