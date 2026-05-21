%% GLORYS annual means and isopycnal interpolation:
clear all; close all; clc

%% Annual means:

anno        = 1993:2015;

file_grid   ='...\OCCITENS\ORCA025.L75-OCCITENS.003\GRID\1993\ORCA025.L75-OCCITENS.003_y1993m01.1m_gridW.nc';
X           = ncread(file_grid,'nav_lon'); 
Y           = ncread(file_grid,'nav_lat');
Z           = ncread(file_grid,'depthw');
dimlon      = size(X,1); dimlat = size(Y,2); dimdep = length(Z);

w_glorys = zeros(dimlon,dimlat,dimdep,length(anno));

file_path = '...\GLORYS12v1\GLORYSORCA025_w\';
ff = dir(file_path);
for ij = 3:length(ff)-3
    disp(ff(ij).name)
    file = strcat(file_path,ff(ij).name);
    load(file)
    w_glorys(:,:,:,ij-2) = mean(wORCA025_75_m,4);
end

file_out = '...\GLORYS12v1\glorys_glob_025_annual.nc';
create_ncfile_verlev(file_out, X, Y, Z, anno, w_glorys, 'w_glorys', 'glorys12v1 (ORCA025 horizontal grid interpolated) vertical velocities');


%% Isopycnal interpolation mean neutral density --------------------------------------------
anno        = 1993:2015;
% Mean density:
file_grid   ='...\GLORYS12v1\GLORYSORCA025_sigma\mercatorglorys4_annual_GLOB_sigma_1993.nc';
X           = ncread(file_grid,'nav_lon'); 
Y           = ncread(file_grid,'nav_lat');
Z           = ncread(file_grid,'deptht');
dimlon      = size(X,1); dimlat = size(Y,2); dimdep = length(Z);

sigma_glorys = zeros(dimlon,dimlat,dimdep,length(anno));

file_path = '...\GLORYS12v1\GLORYSORCA025_sigma\';
ff = dir(file_path);
for ij = 3:length(ff)
    disp(ff(ij).name)
    file = strcat(file_path,ff(ij).name);
    sigma_glorys(:,:,:,ij-2) = squeeze(ncread(file,'vosigntr'));
end
sigmam = mean(sigma_glorys,4);

% Create .nc file ----------------------------------------------------
file_out = '...\GLORYS12v1\glorys_glob_025_annual_sigmam.nc';
nccreate(file_out,'longitude',...
         'Dimensions', {'x',dimlon,'y',dimlat},'Format','netcdf4','Datatype','single')
nccreate(file_out,'latitude',...
         'Dimensions', {'x',dimlon,'y',dimlat},'Format','netcdf4','Datatype','single')
nccreate(file_out,'depth',...
         'Dimensions', {'lev',length(Z)},'Format','netcdf4','Datatype','single')
nccreate(file_out,'signtr',...
         'Dimensions', {'x',dimlon,'y',dimlat,'lev',length(Z)},'Format','netcdf4','Datatype','double')

% Write variables
ncwrite(file_out,'longitude',X)
ncwrite(file_out,'latitude',Y)    
ncwrite(file_out,'depth',squeeze(Z)) 
ncwrite(file_out,'signtr',sigmam) 

% Global Attributes
ncwriteatt(file_out,'/','description',                  'GLORYS netural density mean');
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
ncwriteatt(file_out,'depth',    'long_name',    'Isopycnal levels')
ncwriteatt(file_out,'depth',    'units',        'kg m^-3')
ncwriteatt(file_out,'depth',    'positive',     'down')

%quasi-geostrophic velocities
ncwriteatt(file_out,'signtr',  'axis',         'X,Y,X,T')
ncwriteatt(file_out,'signtr',  'long_name',    'neutral density')
ncwriteatt(file_out,'signtr',  'units','kg/m**3')

%%
disp('Isopycnal interpolation:')

% Isopycnal definition ------------------------------------------------
ind         = [1:54];
isop        = log(ind)./1.26+25;
isopl1      = [21:0.25:26];
isop_uni    = [isopl1 isop(5:end)];
isopl       = isop_uni;
% ----------------------------------------------------------------------

filo = '...\GLORYS12v1\glorys_glob_025_annual.nc';
count = [Inf,Inf,Inf,1];

% Neutral density loading:

filesn  = '...\GLORYS12v1\glorys_glob_025_annual_sigmam.nc';
sigma   = ncread(filesn,'signtr');

for an = 1:length(anno)
    disp(['Year: ',num2str(anno(an))])
    
    start = [1,1,1,an];

    % w loading:

    w_an    = ncread(filo,'w_glorys',start,count);

    [w_isop(:,:,:,an), h_isop(:,:,:,an)] = isop_interp(w_an,dimlon,dimlat,dimdep,isopl,Z,sigma);
end

%% Create .nc file ----------------------------------------------------
file_out = '...\GLORYS12v1\glorys_glob_025_annual_isolevm.nc';
nccreate(file_out,'longitude',...
         'Dimensions', {'x',dimlon,'y',dimlat},'Format','netcdf4','Datatype','single')
nccreate(file_out,'latitude',...
         'Dimensions', {'x',dimlon,'y',dimlat},'Format','netcdf4','Datatype','single')
nccreate(file_out,'isolev',...
         'Dimensions', {'lev',length(isopl)},'Format','netcdf4','Datatype','single')
nccreate(file_out,'time',...
         'Dimensions', {'time',length(anno)},'Format','netcdf4','Datatype','single')
nccreate(file_out,'w_glorys_isolev',...
         'Dimensions', {'x',dimlon,'y',dimlat,'lev',length(isopl),'time',length(anno)},'Format','netcdf4','Datatype','double')
nccreate(file_out,'h_isolev',...
         'Dimensions', {'x',dimlon,'y',dimlat,'lev',length(isopl),'time',length(anno)},'Format','netcdf4','Datatype','double')

% Write variables
ncwrite(file_out,'longitude',X)
ncwrite(file_out,'latitude',Y)    
ncwrite(file_out,'isolev',squeeze(isopl)) 
ncwrite(file_out,'time',anno') 
ncwrite(file_out,'w_glorys_isolev',w_isop) 

% Global Attributes
ncwriteatt(file_out,'/','description',                  'GLORYS vertical velocities');
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

%quasi-geostrophic velocities
ncwriteatt(file_out,'w_glorys_isolev',  'axis',         'X,Y,X,T')
ncwriteatt(file_out,'w_glorys_isolev',  'long_name',    'total vertical velocities')
ncwriteatt(file_out,'w_glorys_isolev',  'units','m/s')

%depth isopycnal surfaces
ncwriteatt(file_out,'h_isolev',  'axis',         'X,Y,X,T')
ncwriteatt(file_out,'h_isolev',  'long_name',    'depth isopycnal surfaces')
ncwriteatt(file_out,'h_isolev',  'units','m')
