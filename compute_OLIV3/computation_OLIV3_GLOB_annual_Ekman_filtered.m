% OLIV3 computation

clear all; close all; clc

% GRID -----------------------------------------------------------------
P           = 'E:\ARMOR3D\1993\dataset-armor-3d-rep-weekly_19930106T1200Z_P20201023T1646Z.nc';
LON         = ncread(P,'longitude');
LAT         = ncread(P,'latitude');
PROF        = ncread(P,'depth');
%LONc        = circshift(LON,length(LON)/2,1);
[LATm,LONm] = meshgrid(double(LAT),double(LON));
dif_zt      = single(abs(PROF(1:end-1)-PROF(2:end)));                       % vertical scale factor

% Parameters -----------------------------------------------------------
rho0        = 1025;
Rt          = 6371229;
anno        = 1993:2019;
fcO         = 1.454441043328608e-04;                                        % 2*rotation_rate
f           = fcO.*(sin(LATm*pi/180)); 

%% ERA5 annual wind stress means ----------------------------------------
file_ERA5   = 'C:\Users\yago_\Documents\LOCEAN\Data\ERA5\ERA5_neutral_wind_and_mean_wind_stress_1993_2020_6_05_22.nc';
lon         = ncread(file_ERA5,'longitude');
lat         = ncread(file_ERA5,'latitude');
time        = ncread(file_ERA5,'time');
taux_ERA    = ncread(file_ERA5,'metss');
tauy_ERA    = ncread(file_ERA5,'mntss');

tauy_ERAr   = reshape(tauy_ERA,[size(tauy_ERA,1) size(tauy_ERA,2) 12 28]); 
taux_ERAr   = reshape(taux_ERA,[size(taux_ERA,1) size(taux_ERA,2) 12 28]);
tauxm       = squeeze(mean(taux_ERAr,3)); 
tauym       = squeeze(mean(tauy_ERAr,3)); 

clear taux_ERAr tauy_ERAr taux_ERA tauy_ERA

%lonc            = lon;                                                      % Circshift
%lonc(lon> 180)  = lon-360;
%lonc            = circshift(lonc,719);

[latm,lonm]     = meshgrid(lat,lon);
fe              = fcO.*(sin(latm*pi/180)); 

% Annual w_GLVB computation --------------------------------------------
for an = 1:length(anno)
    disp(['Year: ',num2str(anno(an))])

    % Geostrophic data loading:
    filev = sprintf('C:\\Users\\yago_\\Documents\\LOCEAN\\Data\\ARMOR3D\\annual_means\\dataset-armor-3d-rep-yearly_%04d.nc',anno(an));
    vgo     = ncread(filev,'vg');
    
    % betav computation:
    df      = f(:,2:end)-f(:,1:end-1);
    dy      = abs(LATm(:,2:end)-LATm(:,1:end-1))*Rt*pi/180;
    
    vgoyc   = circshift(vgo,length(LON)/2,1);
    v       = vgo(:,:,1:end);                       %clear vgoyc
    vj      = 0.5*(v(:,2:end,:)+v(:,1:end-1,:));        clear v
    dfdy    = df./dy;                                   clear df dy
    
    bb      = (vj.*dfdy);                             clear vj dfdy3D   

    % Ekman pumping computation:
    tauymc          = tauym(:,:,an); 
    tauxmc          = tauxm(:,:,an);

    dtauy           = tauymc(2:end,:)./fe(2:end,:)-tauymc(1:end-1,:)./fe(1:end-1,:);        % Zonal gradient tau_y
    dx              = -Rt.*abs(lonm(2:end,:)-lonm(1:end-1,:)).*pi/180;
    dtauydx_f       = dtauy./dx;
      
    dtaux           = (tauxmc(:,2:end)./fe(:,2:end))-(tauxmc(:,1:end-1)./fe(:,1:end-1));    % Meridional gradient tau_x
    dy              = -Rt*abs(latm(:,2:end)-latm(:,1:end-1))*pi/180;
    dtauxdy_f       = dtaux(:,:)./dy(:,:);

    w_Ek            = -(1/rho0)*(dtauydx_f(:,1:end-1)+dtauxdy_f(1:end-1,:));

    % Interpolation ERA5 -> ARMOR3D
    [lonm1,latm1]   = meshgrid(double(lon(1:end-1)),double(lat(1:end-1)));
    w_Ek            = double(w_Ek);
    
    w_Ekint         = interp2(lonm1,latm1,w_Ek',LONm,LATm);
    w_Ekint_masked  = w_Ekint; w_Ekint_masked(isnan(vgo(:,:,1))==1) = NaN;

    %w_ref           = w_Ekint_masked;

    % Filtered Ekman pumping
    w_Ekint_masked(w_Ekint_masked < -10^-5) = 0;
    nfilt           = 4;
    w_Ekf           = smooth2a(double(squeeze(w_Ekint_masked)),nfilt);
    w_ref           = w_Ekf;

    % Vertical integration:
    ff      = (f(:,2:end)+f(:,1:end-1))/2;
    dif_lvb = (1./ff).*(bb(:,:,1:end-1)+bb(:,:,2:end))/2;

    sumlvb          = zeros(size(LON,1),size(LAT,1),length(PROF));
    sumlvb(:,:,1)   = w_ref;

    for i = 1:size(LON,1)-1
        for j = 1:size(LAT,1)-1
            for k = 1:length(PROF)-1
                sumlvb(i,j,k+1) = sumlvb(i,j,k)-dif_lvb(i,j,k)*dif_zt(k);
            end
        end
    end

    % Annual w_GLVB:
    ww_t(:,:,:,an)  = sumlvb;
    
end

% Create .nc file ----------------------------------------------------
file_out = 'C:\\Users\\yago_\\Documents\\LOCEAN\\Data\\OLIV3\\oliv3_glob_025_ekf_annual_verlev.nc';
nccreate(file_out,'longitude',...
         'Dimensions', {'x',size(LONm,1),'y',size(LATm,2)},'Format','netcdf4','Datatype','single')
nccreate(file_out,'latitude',...
         'Dimensions', {'x',size(LONm,1),'y',size(LATm,2)},'Format','netcdf4','Datatype','single')
nccreate(file_out,'depth',...
         'Dimensions', {'deptht',size(PROF,1)},'Format','netcdf4','Datatype','single')
nccreate(file_out,'time',...
         'Dimensions', {'time',length(anno)},'Format','netcdf4','Datatype','single')
nccreate(file_out,'w_oliv3',...
         'Dimensions', {'x',size(LONm,1),'y',size(LATm,2),'deptht',size(PROF,1),'time',length(anno)},'Format','netcdf4','Datatype','double')

% Write variables
ncwrite(file_out,'longitude',LONm)
ncwrite(file_out,'latitude',LATm)    
ncwrite(file_out,'depth',PROF) 
ncwrite(file_out,'time',anno') 
ncwrite(file_out,'w_oliv3',ww_t) 

% Global Attributes
ncwriteatt(file_out,'/','description',                  'Beta-plane geostrophic vertical velocities');
ncwriteatt(file_out,'/','Boundary-condition input data','ERA5 wind stress 2deg filtered');
ncwriteatt(file_out,'/','v_g input data',               'ARMOR3D meridional geostrophic velocities');
ncwriteatt(file_out,'/','output_frequency',             '1yr');
    
%longitude
ncwriteatt(file_out,'longitude','axis',         'X')
ncwriteatt(file_out,'longitude','standard_name','longitude')
ncwriteatt(file_out,'longitude','long_name',    'Longitude')
ncwriteatt(file_out,'longitude','units',        'degrees_east')

%latitude
ncwriteatt(file_out,'latitude', 'axis',         'Y')
ncwriteatt(file_out,'latitude', 'standard_name','latitude')
ncwriteatt(file_out,'latitude', 'long_name',    'Latitude')
ncwriteatt(file_out,'latitude', 'units',        'degrees_north')
    
%deptht
ncwriteatt(file_out,'depth',    'axis',         'Z')
ncwriteatt(file_out,'depth',    'long_name',    'Vertical levels')
ncwriteatt(file_out,'depth',    'units',        'm')
ncwriteatt(file_out,'depth',    'positive',     'down')

%time
ncwriteatt(file_out,'time',     'axis',         'T')
ncwriteatt(file_out,'time',     'long_name',    'year')
ncwriteatt(file_out,'time',     'units',        'year')

%beta-plane geostrophic velocities
ncwriteatt(file_out,'w_oliv3',  'axis',         'X,Y,X,T')
ncwriteatt(file_out,'w_oliv3',  'long_name',    'beta-plane geostrophic vertical velocities')
ncwriteatt(file_out,'w_oliv3',  'units','m/s')

%% Isopycnal interpolation yearly neutral density --------------------------------------------

disp('Isopycnal interpolation:')
anno        = 1993:2019;

% Isopycnal definition ------------------------------------------------
ind         = [1:54];
isop        = log(ind)./1.26+25;
isopl1      = [21:0.25:26];
isop_uni    = [isopl1 isop(5:end)];
isopl       = isop_uni;
% ----------------------------------------------------------------------

for an = 1:length(anno)
    disp(['Year: ',num2str(anno(an))])

    start = [1,1,1,an];
    count = [Inf,Inf,Inf,1];

    % OLIV3 loading:

    filo    = 'C:\\Users\\yago_\\Documents\\LOCEAN\\Data\\OLIV3\\oliv3_glob_025_ekf_annual_verlev.nc';
    w_an    = ncread(filo,'w_oliv3',start,count);

    % Neutral density loading:
    filesn  = sprintf('C:\\Users\\yago_\\Documents\\LOCEAN\\Data\\ARMOR3D\\neutral_density\\dataset-armor-3d-rep-yearly-signtr_%04d.nc',anno(an));
    sigma   = ncread(filesn,'signtr');
    
    dimlon = size(w_an,1); dimlat = size(w_an,2); dimdep = size(w_an,3);

    [w_isop(:,:,:,an), h_isop(:,:,:,an)] = isop_interp(w_an,dimlon,dimlat,dimdep,isopl,PROF,sigma);
end

% Create .nc file ----------------------------------------------------
file_out = 'C:\\Users\\yago_\\Documents\\LOCEAN\\Data\\OLIV3\\oliv3_glob_025_ekf_annual_isolevy.nc';
nccreate(file_out,'longitude',...
         'Dimensions', {'x',size(LONm,1),'y',size(LATm,2)},'Format','netcdf4','Datatype','single')
nccreate(file_out,'latitude',...
         'Dimensions', {'x',size(LONm,1),'y',size(LATm,2)},'Format','netcdf4','Datatype','single')
nccreate(file_out,'isolev',...
         'Dimensions', {'lev',length(isopl)},'Format','netcdf4','Datatype','single')
nccreate(file_out,'time',...
         'Dimensions', {'time',length(anno)},'Format','netcdf4','Datatype','single')
nccreate(file_out,'w_oliv3_isolev',...
         'Dimensions', {'x',size(LONm,1),'y',size(LATm,2),'lev',length(isopl),'time',length(anno)},'Format','netcdf4','Datatype','double')

% Write variables
ncwrite(file_out,'longitude',LONm)
ncwrite(file_out,'latitude',LATm)    
ncwrite(file_out,'isolev',squeeze(isopl)) 
ncwrite(file_out,'time',anno') 
ncwrite(file_out,'w_oliv3_isolev',w_isop) 

% Global Attributes
ncwriteatt(file_out,'/','description',                  'Beta-plane geostrophic vertical velocities');
ncwriteatt(file_out,'/','Boundary-condition input data','ERA5 wind stress');
ncwriteatt(file_out,'/','v_g input data',               'ARMOR3D meridional geostrophic velocities');
ncwriteatt(file_out,'/','output_frequency',             '1yr');
    
%longitude
ncwriteatt(file_out,'longitude','axis',         'X')
ncwriteatt(file_out,'longitude','standard_name','longitude')
ncwriteatt(file_out,'longitude','long_name',    'Longitude')
ncwriteatt(file_out,'longitude','units',        'degrees_east')

%latitude
ncwriteatt(file_out,'latitude', 'axis',         'Y')
ncwriteatt(file_out,'latitude', 'standard_name','latitude')
ncwriteatt(file_out,'latitude', 'long_name',    'Latitude')
ncwriteatt(file_out,'latitude', 'units',        'degrees_north')
    
%deptht
ncwriteatt(file_out,'isolev',    'axis',         'Z')
ncwriteatt(file_out,'isolev',    'long_name',    'Isopycnal levels')
ncwriteatt(file_out,'isolev',    'units',        'kg m^-3')
ncwriteatt(file_out,'isolev',    'positive',     'down')

%time
ncwriteatt(file_out,'time',     'axis',         'T')
ncwriteatt(file_out,'time',     'long_name',    'year')
ncwriteatt(file_out,'time',     'units',        'year')

%beta-plane geostrophic velocities
ncwriteatt(file_out,'w_oliv3_isolev',  'axis',         'X,Y,X,T')
ncwriteatt(file_out,'w_oliv3_isolev',  'long_name',    'beta-plane geostrophic vertical velocities')
ncwriteatt(file_out,'w_oliv3_isolev',  'units','m/s')


% Isopycnal interpolation 27yr mean neutral density --------------------------------------------

disp('Isopycnal interpolation:')
anno        = 1993:2019;

% Isopycnal definition ------------------------------------------------
ind         = [1:54];
isop        = log(ind)./1.26+25;
isopl1      = [21:0.25:26];
isop_uni    = [isopl1 isop(5:end)];
isopl       = isop_uni;
% ----------------------------------------------------------------------

% Neutral density loading:
filesn  = 'C:\\Users\\yago_\\Documents\\LOCEAN\\Data\\ARMOR3D\\neutral_density\\dataset-armor-3d-rep-27yr-signtr.nc';
sigma   = ncread(filesn,'signtr');

for an = 1:length(anno)
    disp(['Year: ',num2str(anno(an))])

    start = [1,1,1,an];
    count = [Inf,Inf,Inf,1];

    % OLIV3 loading:

    filo    = 'C:\\Users\\yago_\\Documents\\LOCEAN\\Data\\OLIV3\\oliv3_glob_025_ekf_annual_verlev.nc';
    w_an    = ncread(filo,'w_oliv3',start,count);
    
    dimlon = size(w_an,1); dimlat = size(w_an,2); dimdep = size(w_an,3);

    [w_isop(:,:,:,an), h_isop(:,:,:,an)] = isop_interp(w_an,dimlon,dimlat,dimdep,isopl,PROF,sigma);
end

% Create .nc file ----------------------------------------------------
file_out = 'C:\\Users\\yago_\\Documents\\LOCEAN\\Data\\OLIV3\\oliv3_glob_025_ekf_annual_isolevm.nc';
nccreate(file_out,'longitude',...
         'Dimensions', {'x',size(LONm,1),'y',size(LATm,2)},'Format','netcdf4','Datatype','single')
nccreate(file_out,'latitude',...
         'Dimensions', {'x',size(LONm,1),'y',size(LATm,2)},'Format','netcdf4','Datatype','single')
nccreate(file_out,'isolev',...
         'Dimensions', {'lev',length(isopl)},'Format','netcdf4','Datatype','single')
nccreate(file_out,'time',...
         'Dimensions', {'time',length(anno)},'Format','netcdf4','Datatype','single')
nccreate(file_out,'w_oliv3_isolev',...
         'Dimensions', {'x',size(LONm,1),'y',size(LATm,2),'lev',length(isopl),'time',length(anno)},'Format','netcdf4','Datatype','double')
nccreate(file_out,'h_isolev',...
         'Dimensions', {'x',size(LONm,1),'y',size(LATm,2),'lev',length(isopl),'time',length(anno)},'Format','netcdf4','Datatype','double')


% Write variables
ncwrite(file_out,'longitude',LONm)
ncwrite(file_out,'latitude',LATm)    
ncwrite(file_out,'isolev',squeeze(isopl)) 
ncwrite(file_out,'time',anno') 
ncwrite(file_out,'w_oliv3_isolev',w_isop) 
ncwrite(file_out,'h_isolev',h_isop) 

% Global Attributes
ncwriteatt(file_out,'/','description',                  'Beta-plane geostrophic vertical velocities');
ncwriteatt(file_out,'/','Boundary-condition input data','ERA5 wind stress');
ncwriteatt(file_out,'/','v_g input data',               'ARMOR3D meridional geostrophic velocities');
ncwriteatt(file_out,'/','output_frequency',             '1yr');
    
%longitude
ncwriteatt(file_out,'longitude','axis',         'X')
ncwriteatt(file_out,'longitude','standard_name','longitude')
ncwriteatt(file_out,'longitude','long_name',    'Longitude')
ncwriteatt(file_out,'longitude','units',        'degrees_east')

%latitude
ncwriteatt(file_out,'latitude', 'axis',         'Y')
ncwriteatt(file_out,'latitude', 'standard_name','latitude')
ncwriteatt(file_out,'latitude', 'long_name',    'Latitude')
ncwriteatt(file_out,'latitude', 'units',        'degrees_north')
    
%deptht
ncwriteatt(file_out,'isolev',    'axis',         'Z')
ncwriteatt(file_out,'isolev',    'long_name',    'Isopycnal levels')
ncwriteatt(file_out,'isolev',    'units',        'kg m^-3')
ncwriteatt(file_out,'isolev',    'positive',     'down')

%time
ncwriteatt(file_out,'time',     'axis',         'T')
ncwriteatt(file_out,'time',     'long_name',    'year')
ncwriteatt(file_out,'time',     'units',        'year')

%beta-plane geostrophic velocities
ncwriteatt(file_out,'w_oliv3_isolev',  'axis',         'X,Y,X,T')
ncwriteatt(file_out,'w_oliv3_isolev',  'long_name',    'beta-plane geostrophic vertical velocities')
ncwriteatt(file_out,'w_oliv3_isolev',  'units','m/s')

%depth isopycnal surfaces
ncwriteatt(file_out,'h_isolev',  'axis',         'X,Y,X,T')
ncwriteatt(file_out,'h_isolev',  'long_name',    'depth isopycnal surfaces')
ncwriteatt(file_out,'h_isolev',  'units','m')

% Smoothing filter -------------------------------------------

anno        = 1993:2019;

% Mean isopycnal surface:
nfilt       = 5*4/2;

for an = 1:length(anno)
    disp(['Year: ',num2str(anno(an))])

    start = [1,1,1,an];
    count = [Inf,Inf,Inf,1];

    % OLIV3 loading:

    filo    = 'C:\\Users\\yago_\\Documents\\LOCEAN\\Data\\OLIV3\\oliv3_glob_025_ekf_annual_isolevm.nc';
    w_an    = ncread(filo,'w_oliv3_isolev',start,count);
    
    if an == 1
        lat         = ncread(filo,'latitude'); 
        lon         = ncread(filo,'longitude');
        isolev      = ncread(filo,'isolev');
    end
    
    dimlon = size(w_an,1); dimlat = size(w_an,2); dimdep = size(w_an,3);
    
    for k = 1:dimdep
        wk              = squeeze(w_an(:,:,k));
        wf(:,:,k,an)    = smooth2a(wk,nfilt);
    end
end

% Create .nc file ----------------------------------------------------
file_out = 'C:\\Users\\yago_\\Documents\\LOCEAN\\Data\\OLIV3\\oliv3_glob_025_ekf_annual_isolevm_5filtr.nc';
nccreate(file_out,'longitude',...
         'Dimensions', {'x',size(LONm,1),'y',size(LATm,2)},'Format','netcdf4','Datatype','single')
nccreate(file_out,'latitude',...
         'Dimensions', {'x',size(LONm,1),'y',size(LATm,2)},'Format','netcdf4','Datatype','single')
nccreate(file_out,'isolev',...
         'Dimensions', {'lev',length(isolev)},'Format','netcdf4','Datatype','single')
nccreate(file_out,'time',...
         'Dimensions', {'time',length(anno)},'Format','netcdf4','Datatype','single')
nccreate(file_out,'w_oliv3_isolevf',...
         'Dimensions', {'x',size(LONm,1),'y',size(LATm,2),'lev',length(isolev),'time',length(anno)},'Format','netcdf4','Datatype','double')

% Write variables
ncwrite(file_out,'longitude',LONm)
ncwrite(file_out,'latitude',LATm)    
ncwrite(file_out,'isolev',isolev) 
ncwrite(file_out,'time',anno') 
ncwrite(file_out,'w_oliv3_isolevf',wf) 

% Global Attributes
ncwriteatt(file_out,'/','description',                  'Beta-plane geostrophic vertical velocities');
ncwriteatt(file_out,'/','Boundary-condition input data','ERA5 wind stress');
ncwriteatt(file_out,'/','v_g input data',               'ARMOR3D meridional geostrophic velocities');
ncwriteatt(file_out,'/','Isopycnal interpolation',      'computed using ARMOR3D 27yr means of neutral density');
ncwriteatt(file_out,'/','output_frequency',             '1yr');
    
%longitude
ncwriteatt(file_out,'longitude','axis',         'X')
ncwriteatt(file_out,'longitude','standard_name','longitude')
ncwriteatt(file_out,'longitude','long_name',    'Longitude')
ncwriteatt(file_out,'longitude','units',        'degrees_east')

%latitude
ncwriteatt(file_out,'latitude', 'axis',         'Y')
ncwriteatt(file_out,'latitude', 'standard_name','latitude')
ncwriteatt(file_out,'latitude', 'long_name',    'Latitude')
ncwriteatt(file_out,'latitude', 'units',        'degrees_north')
    
%deptht
ncwriteatt(file_out,'isolev',    'axis',         'Z')
ncwriteatt(file_out,'isolev',    'long_name',    'Isopycnal levels')
ncwriteatt(file_out,'isolev',    'units',        'kg m^-3')
ncwriteatt(file_out,'isolev',    'positive',     'down')

%time
ncwriteatt(file_out,'time',     'axis',         'T')
ncwriteatt(file_out,'time',     'long_name',    'year')
ncwriteatt(file_out,'time',     'units',        'year')

%beta-plane geostrophic velocities
ncwriteatt(file_out,'w_oliv3_isolevf',  'axis',         'X,Y,X,T')
ncwriteatt(file_out,'w_oliv3_isolevf',  'long_name',    'beta-plane geostrophic vertical velocities smoothed with 5deg filter')
ncwriteatt(file_out,'w_oliv3_isolevf',  'units','m/s')


clear wf
% Yearly isopycnal surface:
nfilt       = 5*4/2;

for an = 1:length(anno)
    disp(['Year: ',num2str(anno(an))])

    start = [1,1,1,an];
    count = [Inf,Inf,Inf,1];

    % OLIV3 loading:

    filo    = 'C:\\Users\\yago_\\Documents\\LOCEAN\\Data\\OLIV3\\oliv3_glob_025_ekf_annual_isolevy.nc';
    w_an    = ncread(filo,'w_oliv3_isolev',start,count);
    
    dimlon = size(w_an,1); dimlat = size(w_an,2); dimdep = size(w_an,3);
    
    for k = 1:dimdep
        wk              = squeeze(w_an(:,:,k));
        wf(:,:,k,an)    = smooth2a(wk,nfilt);
    end
end

% Create .nc file ----------------------------------------------------
file_out = 'C:\\Users\\yago_\\Documents\\LOCEAN\\Data\\OLIV3\\oliv3_glob_025_ekf_annual_isolevy_5filtr.nc';
nccreate(file_out,'longitude',...
         'Dimensions', {'x',size(LONm,1),'y',size(LATm,2)},'Format','netcdf4','Datatype','single')
nccreate(file_out,'latitude',...
         'Dimensions', {'x',size(LONm,1),'y',size(LATm,2)},'Format','netcdf4','Datatype','single')
nccreate(file_out,'isolev',...
         'Dimensions', {'lev',length(isolev)},'Format','netcdf4','Datatype','single')
nccreate(file_out,'time',...
         'Dimensions', {'time',length(anno)},'Format','netcdf4','Datatype','single')
nccreate(file_out,'w_oliv3_isolevf',...
         'Dimensions', {'x',size(LONm,1),'y',size(LATm,2),'lev',length(isolev),'time',length(anno)},'Format','netcdf4','Datatype','double')

% Write variables
ncwrite(file_out,'longitude',LONm)
ncwrite(file_out,'latitude',LATm)    
ncwrite(file_out,'isolev',isolev) 
ncwrite(file_out,'time',anno') 
ncwrite(file_out,'w_oliv3_isolevf',wf) 

% Global Attributes
ncwriteatt(file_out,'/','description',                  'Beta-plane geostrophic vertical velocities');
ncwriteatt(file_out,'/','Boundary-condition input data','ERA5 wind stress');
ncwriteatt(file_out,'/','v_g input data',               'ARMOR3D meridional geostrophic velocities');
ncwriteatt(file_out,'/','Isopycnal interpolation',      'computed using ARMOR3D annual means of neutral density');
ncwriteatt(file_out,'/','output_frequency',             '1yr');
    
%longitude
ncwriteatt(file_out,'longitude','axis',         'X')
ncwriteatt(file_out,'longitude','standard_name','longitude')
ncwriteatt(file_out,'longitude','long_name',    'Longitude')
ncwriteatt(file_out,'longitude','units',        'degrees_east')

%latitude
ncwriteatt(file_out,'latitude', 'axis',         'Y')
ncwriteatt(file_out,'latitude', 'standard_name','latitude')
ncwriteatt(file_out,'latitude', 'long_name',    'Latitude')
ncwriteatt(file_out,'latitude', 'units',        'degrees_north')
    
%deptht
ncwriteatt(file_out,'isolev',    'axis',         'Z')
ncwriteatt(file_out,'isolev',    'long_name',    'Isopycnal levels')
ncwriteatt(file_out,'isolev',    'units',        'kg m^-3')
ncwriteatt(file_out,'isolev',    'positive',     'down')

%time
ncwriteatt(file_out,'time',     'axis',         'T')
ncwriteatt(file_out,'time',     'long_name',    'year')
ncwriteatt(file_out,'time',     'units',        'year')

%beta-plane geostrophic velocities
ncwriteatt(file_out,'w_oliv3_isolevf',  'axis',         'X,Y,X,T')
ncwriteatt(file_out,'w_oliv3_isolevf',  'long_name',    'beta-plane geostrophic vertical velocities smoothed with 5deg filter')
ncwriteatt(file_out,'w_oliv3_isolevf',  'units','m/s')


%% FIGURE test

cont2 = [-1.8 -1.6 -1.4 -1.2 -1 -0.8 -0.6 -0.4 -0.2 0 0.2 0.4 0.6 0.8 1 1.2]*10^(-5);
cont3 = [-20 -1 -0.8 -0.6 -0.4 -0.2 0 0.2 0.4 0.6 0.8 1]*10^(-5);
cont10 = [-1 -0.8 -0.6 -0.4 -0.2 0 0.2 0.4 0.6 0.8 1]*10^(-6);

RB1= cbrewer('div', 'RdBu',length(-20:1:20)); 

RB1=RB1'; RB1=fliplr(RB1);RB1=RB1';

RB1_sel = [RB1(1:2:9,:);RB1(33:2:41,:)];
RB11c = [RB1(10:2:18,:);RB1(24:2:32,:)];

% Mapa de color variable

RBens = zeros(100,3);

for i = 1:3
RBens(1:10,i) = RB1_sel(1,i);
RBens(11:20,i) = RB1_sel(2,i);
RBens(21:30,i) = RB1_sel(3,i);
RBens(31:40,i) = RB1_sel(4,i);
RBens(41:45,i) = RB1_sel(5,i);
RBens(46:55,i) = RB11c(:,i);
RBens(56:60,i) = RB1_sel(6,i);
RBens(61:70,i) = RB1_sel(7,i);
RBens(71:80,i) = RB1_sel(8,i);
RBens(81:90,i) = RB1_sel(9,i);
RBens(91:100,i) = RB1_sel(10,i);
end

contens = [cont3(1:6) cont10 cont3(8:end)];
contens2 = [cont3(2:6) cont10 cont3(8:end)];
%%
for k2 = 21%[19 21 24 50 57 63 69]

    figure
    subplot(2,2,1)
    colormap(RBens)

    % Yearly isopycnals no filter;
    filo    = 'C:\\Users\\yago_\\Documents\\LOCEAN\\Data\\OLIV3\\oliv3_glob_025_ekf_annual_isolevy.nc';
    wm    = mean(ncread(filo,'w_oliv3_isolev'),4);
    %wm      = mean(w_an,4);

    [c1 h1] = contourf(squeeze(wm(:,:,k2))',contens,'LineColor','none');

    h1 = colorbar('southoutside');
    set(h1,'xtick',cont2); 
    set(get(h1,'xlabel'),'string','m/s');
    caxis([-10 10]*10^-6)

    title(['Yearly isopycnals no filter \sigma = ',num2str(round(isopl(k2),2)),' kg/m^{3}']);
    
    subplot(2,2,2)
    colormap(RBens)
    % Mean isopycnals no filter;
    filo    = 'C:\\Users\\yago_\\Documents\\LOCEAN\\Data\\OLIV3\\oliv3_glob_025_ekf_annual_isolevm.nc';
    wm    = mean(ncread(filo,'w_oliv3_isolev'),4);
    %wm      = mean(w_an,4);

    [c1 h1] = contourf(squeeze(wm(:,:,k2))',contens,'LineColor','none');
    
    h1 = colorbar('southoutside');
    set(h1,'xtick',cont2); 
    set(get(h1,'xlabel'),'string','m/s');
    caxis([-10 10]*10^-6)

    title(['Mean isopycnals no filter \sigma = ',num2str(round(isopl(k2),2)),' kg/m^{3}']);
    
    subplot(2,2,3)
    colormap(RBens)
    % Yearly isopycnals filter;
    filo    = 'C:\\Users\\yago_\\Documents\\LOCEAN\\Data\\OLIV3\\oliv3_glob_025_ekf_annual_isolevy_5filtr.nc';
    wm    = mean(ncread(filo,'w_oliv3_isolevf'),4);
    %wm      = mean(w_an,4);

    [c1 h1] = contourf(squeeze(wm(:,:,k2))',contens,'LineColor','none');

    h1 = colorbar('southoutside');
    set(h1,'xtick',cont2); 
    set(get(h1,'xlabel'),'string','m/s');
    caxis([-10 10]*10^-6)

    title(['Yearly isopycnals filter \sigma = ',num2str(round(isopl(k2),2)),' kg/m^{3}']);
    
    subplot(2,2,4)
    colormap(RBens)
    % Mean isopycnals filter;
    filo    = 'C:\\Users\\yago_\\Documents\\LOCEAN\\Data\\OLIV3\\oliv3_glob_025_ekf_annual_isolevm_5filtr.nc';
    wm    = mean(ncread(filo,'w_oliv3_isolevf'),4);
    %wm      = mean(w_an,4);

    [c1 h1] = contourf(squeeze(wm(:,:,k2))',contens,'LineColor','none');

    h1 = colorbar('southoutside');
    set(h1,'xtick',cont2); 
    set(get(h1,'xlabel'),'string','m/s');
    caxis([-10 10]*10^-6)

    title(['Mean isopycnals filter \sigma = ',num2str(round(isopl(k2),2)),' kg/m^{3}']);
end

 %% ERA5 annual wind stress means ----------------------------------------
file_ERA5   = 'C:\Users\yago_\Documents\LOCEAN\Data\ERA5\ERA5_neutral_wind_and_mean_wind_stress_1993_2020_6_05_22.nc';
lon         = ncread(file_ERA5,'longitude');
lat         = ncread(file_ERA5,'latitude');
time        = ncread(file_ERA5,'time');
taux_ERA    = ncread(file_ERA5,'metss');
tauy_ERA    = ncread(file_ERA5,'mntss');

tauy_ERAr   = reshape(tauy_ERA,[size(tauy_ERA,1) size(tauy_ERA,2) 12 28]); 
taux_ERAr   = reshape(taux_ERA,[size(taux_ERA,1) size(taux_ERA,2) 12 28]);
tauxm       = squeeze(mean(taux_ERAr,3)); 
tauym       = squeeze(mean(tauy_ERAr,3)); 

clear taux_ERAr tauy_ERAr taux_ERA tauy_ERA

%lonc            = lon;                                                      % Circshift
%lonc(lon> 180)  = lon-360;
%lonc            = circshift(lonc,719);

[latm,lonm]     = meshgrid(lat,lon);
fe              = fcO.*(sin(latm*pi/180)); 

% Annual w_GLVB computation --------------------------------------------
for an = 1:length(anno)
    disp(['Year: ',num2str(anno(an))])
    if an == 1
        % Geostrophic data loading:
        filev = sprintf('C:\\Users\\yago_\\Documents\\LOCEAN\\Data\\ARMOR3D\\annual_means\\dataset-armor-3d-rep-yearly_%04d.nc',anno(an));
        vgo     = ncread(filev,'vg');
    end
    % Ekman pumping computation:
    tauymc          = tauym(:,:,an); 
    tauxmc          = tauxm(:,:,an);

    dtauy           = tauymc(2:end,:)./fe(2:end,:)-tauymc(1:end-1,:)./fe(1:end-1,:);        % Zonal gradient tau_y
    dx              = -Rt.*abs(lonm(2:end,:)-lonm(1:end-1,:)).*pi/180;
    dtauydx_f       = dtauy./dx;
      
    dtaux           = (tauxmc(:,2:end)./fe(:,2:end))-(tauxmc(:,1:end-1)./fe(:,1:end-1));    % Meridional gradient tau_x
    dy              = -Rt*abs(latm(:,2:end)-latm(:,1:end-1))*pi/180;
    dtauxdy_f       = dtaux(:,:)./dy(:,:);

    w_Ek            = -(1/rho0)*(dtauydx_f(:,1:end-1)+dtauxdy_f(1:end-1,:));

    % Interpolation ERA5 -> ARMOR3D
    [lonm1,latm1]   = meshgrid(double(lon(1:end-1)),double(lat(1:end-1)));
    w_Ek            = double(w_Ek);
    
    w_Ekint         = interp2(lonm1,latm1,w_Ek',LONm,LATm);
    w_Ekint_masked  = w_Ekint; w_Ekint_masked(isnan(vgo(:,:,1))==1) = NaN;
    w_Ekint_masked(abs(LATm)<5) = NaN;
    %w_ref           = w_Ekint_masked;

    % Filtered Ekman pumping
    w_Ekint_masked(abs(w_Ekint_masked) > 10^-5) = 0;
    nfilt           = 7;
    w_Ekf           = smooth2a(double(squeeze(w_Ekint_masked)),nfilt);
    w_ref(:,:,an)   = w_Ekf;

end

% Create .nc file ----------------------------------------------------
file_out = 'C:\\Users\\yago_\\Documents\\LOCEAN\\Data\\ERA5\\ekman_pumping_glob_025_filtered_annual.nc';
nccreate(file_out,'longitude',...
         'Dimensions', {'x',size(LONm,1),'y',size(LATm,2)},'Format','netcdf4','Datatype','single')
nccreate(file_out,'latitude',...
         'Dimensions', {'x',size(LONm,1),'y',size(LATm,2)},'Format','netcdf4','Datatype','single')
nccreate(file_out,'time',...
         'Dimensions', {'time',length(anno)},'Format','netcdf4','Datatype','single')
nccreate(file_out,'wek_era5',...
         'Dimensions', {'x',size(LONm,1),'y',size(LATm,2),'time',length(anno)},'Format','netcdf4','Datatype','double')

% Write variables
ncwrite(file_out,'longitude',LONm)
ncwrite(file_out,'latitude',LATm)    
ncwrite(file_out,'time',anno') 
ncwrite(file_out,'wek_era5',w_ref) 

% Global Attributes
ncwriteatt(file_out,'/','description',                  'Ekman pumping vertical velocity');
ncwriteatt(file_out,'/','Boundary-condition input data','ERA5 wind stress 5deg filtered');
ncwriteatt(file_out,'/','output_frequency',             '1yr');
    
%longitude
ncwriteatt(file_out,'longitude','axis',         'X')
ncwriteatt(file_out,'longitude','standard_name','longitude')
ncwriteatt(file_out,'longitude','long_name',    'Longitude')
ncwriteatt(file_out,'longitude','units',        'degrees_east')

%latitude
ncwriteatt(file_out,'latitude', 'axis',         'Y')
ncwriteatt(file_out,'latitude', 'standard_name','latitude')
ncwriteatt(file_out,'latitude', 'long_name',    'Latitude')
ncwriteatt(file_out,'latitude', 'units',        'degrees_north')

%time
ncwriteatt(file_out,'time',     'axis',         'T')
ncwriteatt(file_out,'time',     'long_name',    'year')
ncwriteatt(file_out,'time',     'units',        'year')

%beta-plane geostrophic velocities
ncwriteatt(file_out,'wek_era5',  'axis',         'X,Y,X,T')
ncwriteatt(file_out,'wek_era5',  'long_name',    'Ekman pumping vertical velocities')
ncwriteatt(file_out,'wek_era5',  'units','m/s')

%%
figure
contourf(LONm,LATm,squeeze(mean(w_ref,3)))
