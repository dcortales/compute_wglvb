%% GLORYS annual means and isopycnal interpolation:
clear all; close all; clc

%% Annual means:

anno        = 1993:2015;

file_grid   ='E:\ORCA025.L75-OCCITENS.003\GRID\1993\ORCA025.L75-OCCITENS.003_y1993m01.1m_gridW.nc';
X           = ncread(file_grid,'nav_lon'); 
Y           = ncread(file_grid,'nav_lat');
Z           = ncread(file_grid,'depthw');
dimlon      = size(X,1); dimlat = size(Y,2); dimdep = length(Z);

w_glorys = zeros(dimlon,dimlat,dimdep,length(anno));

file_path = 'E:\GLORYS12v1\GLORYSORCA025_w\';
ff = dir(file_path);
for ij = 3:length(ff)-3
    disp(ff(ij).name)
    file = strcat(file_path,ff(ij).name);
    load(file)
    w_glorys(:,:,:,ij-2) = mean(wORCA025_75_m,4);
end

file_out = 'C:\Users\yago_\Documents\LOCEAN\Data\GLORYS\glorys_glob_025_annual.nc';
create_ncfile_verlev(file_out, X, Y, Z, anno, w_glorys, 'w_glorys', 'glorys12v1 (ORCA025 horizontal grid interpolated) vertical velocities');


%% Isopycnal interpolation mean neutral density --------------------------------------------
anno        = 1993:2015;
% Mean density:
file_grid   ='E:\GLORYS12v1\GLORYSORCA025_sigma\mercatorglorys4_annual_GLOB_sigma_1993.nc';
X           = ncread(file_grid,'nav_lon'); 
Y           = ncread(file_grid,'nav_lat');
Z           = ncread(file_grid,'deptht');
dimlon      = size(X,1); dimlat = size(Y,2); dimdep = length(Z);

sigma_glorys = zeros(dimlon,dimlat,dimdep,length(anno));

file_path = 'E:\GLORYS12v1\GLORYSORCA025_sigma\';
ff = dir(file_path);
for ij = 3:length(ff)
    disp(ff(ij).name)
    file = strcat(file_path,ff(ij).name);
    sigma_glorys(:,:,:,ij-2) = squeeze(ncread(file,'vosigntr'));
end
sigmam = mean(sigma_glorys,4);

% Create .nc file ----------------------------------------------------
file_out = 'C:\Users\yago_\Documents\LOCEAN\Data\GLORYS\glorys_glob_025_annual_sigmam.nc';
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

filo = 'C:\Users\yago_\Documents\LOCEAN\Data\GLORYS\glorys_glob_025_annual.nc';
count = [Inf,Inf,Inf,1];

% Neutral density loading:

filesn  = 'C:\Users\yago_\Documents\LOCEAN\Data\GLORYS\glorys_glob_025_annual_sigmam.nc';
sigma   = ncread(filesn,'signtr');

for an = 1:length(anno)
    disp(['Year: ',num2str(anno(an))])
    
    start = [1,1,1,an];

    % w loading:

    w_an    = ncread(filo,'w_glorys',start,count);

    [w_isop(:,:,:,an), h_isop(:,:,:,an)] = isop_interp(w_an,dimlon,dimlat,dimdep,isopl,Z,sigma);
end

%% Create .nc file ----------------------------------------------------
file_out = 'C:\Users\yago_\Documents\LOCEAN\Data\GLORYS\glorys_glob_025_annual_isolevm.nc';
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


%% Isopycnal interpolation 27yr mean neutral density --------------------------------------------

disp('Isopycnal interpolation:')
anno        = 1993:2019;

% Isopycnal definition ------------------------------------------------
ind         = [1:54];
isop        = log(ind)./1.26+25;
isopl1      = [21:0.25:26];
isop_uni    = [isopl1 isop(5:end)];
isopl       = isop_uni;
% ----------------------------------------------------------------------

file_path   ='E:\ORCA025.L75-OCCITENS.003\GRID\1993\ORCA025.L75-OCCITENS.003_y1993m01.1m_gridW.nc';
X           = ncread(file_path,'nav_lon'); 
Y           = ncread(file_path,'nav_lat');
Z           = ncread(file_path,'depthw');
dimlon      = size(X,1); dimlat = size(Y,2); dimdep = length(Z);

% w occitens

filo    = sprintf('C:\\Users\\yago_\\Documents\\LOCEAN\\Data\\OCCITENS\\w_annual_%04d.nc',anno(1));
w_an    = ncread(filo,'w');

% Neutral density loading:
filesn  = 'C:\\Users\\yago_\\Documents\\LOCEAN\\Data\\OCCITENS\\neutral_density\\ORCA025.L75-OCCITENS.003-signtr_27yr.nc';
sigma   = ncread(filesn,'signtr');

for an = 1:length(anno)
    disp(['Year: ',num2str(anno(an))])

    filo    = sprintf('C:\\Users\\yago_\\Documents\\LOCEAN\\Data\\OCCITENS\\w_annual_%04d.nc',anno(an));
    w_an    = permute(ncread(filo,'w'),[3 2 1]);

    [w_isop(:,:,:,an), h_isop(:,:,:,an)] = isop_interp(w_an,dimlon,dimlat,dimdep,isopl,Z,sigma);
end

% Create .nc file ----------------------------------------------------
file_out = 'C:\\Users\\yago_\\Documents\\LOCEAN\\Data\\OCCITENS\\occitens_glob_025_annual_isolevm.nc';
nccreate(file_out,'longitude',...
         'Dimensions', {'x',dimlon,'y',dimlat},'Format','netcdf4','Datatype','single')
nccreate(file_out,'latitude',...
         'Dimensions', {'x',dimlon,'y',dimlat},'Format','netcdf4','Datatype','single')
nccreate(file_out,'isolev',...
         'Dimensions', {'lev',length(isopl)},'Format','netcdf4','Datatype','single')
nccreate(file_out,'time',...
         'Dimensions', {'time',length(anno)},'Format','netcdf4','Datatype','single')
nccreate(file_out,'w_occitens_isolev',...
         'Dimensions', {'x',dimlon,'y',dimlat,'lev',length(isopl),'time',length(anno)},'Format','netcdf4','Datatype','double')
nccreate(file_out,'h_isolev',...
         'Dimensions', {'x',dimlon,'y',dimlat,'lev',length(isopl),'time',length(anno)},'Format','netcdf4','Datatype','double')

% Write variables
ncwrite(file_out,'longitude',X)
ncwrite(file_out,'latitude',Y)    
ncwrite(file_out,'isolev',squeeze(isopl)) 
ncwrite(file_out,'time',anno') 
ncwrite(file_out,'w_occitens_isolev',w_isop) 

% Global Attributes
ncwriteatt(file_out,'/','description',                  'OCCITENS vertical velocities');
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
ncwriteatt(file_out,'w_occitens_isolev',  'axis',         'X,Y,X,T')
ncwriteatt(file_out,'w_occitens_isolev',  'long_name',    'total vertical velocities')
ncwriteatt(file_out,'w_occitens_isolev',  'units','m/s')

%depth isopycnal surfaces
ncwriteatt(file_out,'h_isolev',  'axis',         'X,Y,X,T')
ncwriteatt(file_out,'h_isolev',  'long_name',    'depth isopycnal surfaces')
ncwriteatt(file_out,'h_isolev',  'units','m')


%% --------------------------------------------------------------------
% OCCITENS w_glvb 
%% --------------------------------------------------------------------

% Isopycnal interpolation yearly neutral density --------------------------------------------

disp('Isopycnal interpolation:')
anno        = 1993:2015;

file_path   ='E:\ORCA025.L75-OCCITENS.003\GRID\1993\ORCA025.L75-OCCITENS.003_y1993m01.1m_gridW.nc';
X           = ncread(file_path,'nav_lon'); 
Y           = ncread(file_path,'nav_lat');
Z           = ncread(file_path,'depthw');
dimlon      = size(X,1); dimlat = size(Y,2); dimdep = length(Z);

% Isopycnal definition ------------------------------------------------
ind         = [1:54];
isop        = log(ind)./1.26+25;
isopl1      = [21:0.25:26];
isop_uni    = [isopl1 isop(5:end)];
isopl       = isop_uni;
% ----------------------------------------------------------------------

for an = 1:length(anno)
    disp(['Year: ',num2str(anno(an))])

    filo    = sprintf('C:\\Users\\yago_\\Documents\\LOCEAN\\Data\\OCCITENS\\wglvb_annual_%04d.nc',anno(an));
    w_an    = permute(ncread(filo,'div'),[3 2 1]);

    % Neutral density loading:
    filesn  = sprintf('C:\\Users\\yago_\\Documents\\LOCEAN\\Data\\OCCITENS\\neutral_density\\ORCA025.L75-OCCITENS.003-signtr_%04d.nc',anno(an));
    sigma   = ncread(filesn,'signtr');

    [w_isop(:,:,:,an), h_isop(:,:,:,an)] = isop_interp(w_an,dimlon,dimlat,dimdep,isopl,Z,sigma);
end

% Create .nc file ----------------------------------------------------
file_out = 'C:\\Users\\yago_\\Documents\\LOCEAN\\Data\\OCCITENS\\occitens_wglvb_glob_025_annual_isolevy.nc';
nccreate(file_out,'longitude',...
         'Dimensions', {'x',dimlon,'y',dimlat},'Format','netcdf4','Datatype','single')
nccreate(file_out,'latitude',...
         'Dimensions', {'x',dimlon,'y',dimlat},'Format','netcdf4','Datatype','single')
nccreate(file_out,'isolev',...
         'Dimensions', {'lev',length(isopl)},'Format','netcdf4','Datatype','single')
nccreate(file_out,'time',...
         'Dimensions', {'time',length(anno)},'Format','netcdf4','Datatype','single')
nccreate(file_out,'wglvb_occitens_isolev',...
         'Dimensions', {'x',dimlon,'y',dimlat,'lev',length(isopl),'time',length(anno)},'Format','netcdf4','Datatype','double')
nccreate(file_out,'h_isolev',...
         'Dimensions', {'x',dimlon,'y',dimlat,'lev',length(isopl),'time',length(anno)},'Format','netcdf4','Datatype','double')

% Write variables
ncwrite(file_out,'longitude',X)
ncwrite(file_out,'latitude',Y)    
ncwrite(file_out,'isolev',squeeze(isopl)) 
ncwrite(file_out,'time',anno') 
ncwrite(file_out,'wglvb_occitens_isolev',w_isop) 

% Global Attributes
ncwriteatt(file_out,'/','description',                  'OCCITENS beta-plane geostrophic vertical velocities');
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
ncwriteatt(file_out,'wglvb_occitens_isolev',  'axis',         'X,Y,X,T')
ncwriteatt(file_out,'wglvb_occitens_isolev',  'long_name',    'beta-plane geostrophic vertical velocities')
ncwriteatt(file_out,'wglvb_occitens_isolev',  'units','m/s')

%depth isopycnal surfaces
ncwriteatt(file_out,'h_isolev',  'axis',         'X,Y,X,T')
ncwriteatt(file_out,'h_isolev',  'long_name',    'depth isopycnal surfaces')
ncwriteatt(file_out,'h_isolev',  'units','m')


%% Isopycnal interpolation 27yr mean neutral density --------------------------------------------

disp('Isopycnal interpolation:')
anno        = 1993:2015;

% Isopycnal definition ------------------------------------------------
ind         = [1:54];
isop        = log(ind)./1.26+25;
isopl1      = [21:0.25:26];
isop_uni    = [isopl1 isop(5:end)];
isopl       = isop_uni;
% ----------------------------------------------------------------------

file_path   ='E:\ORCA025.L75-OCCITENS.003\GRID\1993\ORCA025.L75-OCCITENS.003_y1993m01.1m_gridW.nc';
X           = ncread(file_path,'nav_lon'); 
Y           = ncread(file_path,'nav_lat');
Z           = ncread(file_path,'depthw');
dimlon      = size(X,1); dimlat = size(Y,2); dimdep = length(Z);

% w occitens

filo    = sprintf('C:\\Users\\yago_\\Documents\\LOCEAN\\Data\\OCCITENS\\wglvb_annual_%04d.nc',anno(1));
w_an    = permute(ncread(filo,'div'),[3 2 1]);

% Neutral density loading:
filesn  = 'C:\\Users\\yago_\\Documents\\LOCEAN\\Data\\OCCITENS\\neutral_density\\ORCA025.L75-OCCITENS.003-signtr_27yr.nc';
sigma   = ncread(filesn,'signtr');

for an = 1:length(anno)
    disp(['Year: ',num2str(anno(an))])

    filo    = sprintf('C:\\Users\\yago_\\Documents\\LOCEAN\\Data\\OCCITENS\\wglvb_annual_%04d.nc',anno(an));
    w_an    = permute(ncread(filo,'div'),[3 2 1]);

    [w_isop(:,:,:,an), h_isop] = isop_interp(w_an,dimlon,dimlat,dimdep,isopl,Z,sigma);
end

%% Create .nc file ----------------------------------------------------
file_out = 'E:\OLIV3_figures\OCCITENS\occitens_wglvb_glob_025_annual_isolevm.nc';
nccreate(file_out,'longitude',...
         'Dimensions', {'x',dimlon,'y',dimlat},'Format','netcdf4','Datatype','single')
nccreate(file_out,'latitude',...
         'Dimensions', {'x',dimlon,'y',dimlat},'Format','netcdf4','Datatype','single')
nccreate(file_out,'isolev',...
         'Dimensions', {'lev',length(isopl)},'Format','netcdf4','Datatype','single')
nccreate(file_out,'time',...
         'Dimensions', {'time',length(anno)},'Format','netcdf4','Datatype','single')
nccreate(file_out,'wglvb_occitens_isolev',...
         'Dimensions', {'x',dimlon,'y',dimlat,'lev',length(isopl),'time',length(anno)},'Format','netcdf4','Datatype','double')
nccreate(file_out,'h_isolev',...
         'Dimensions', {'x',dimlon,'y',dimlat,'lev',length(isopl),'time',length(anno)},'Format','netcdf4','Datatype','double')

% Write variables
ncwrite(file_out,'longitude',X)
ncwrite(file_out,'latitude',Y)    
ncwrite(file_out,'isolev',squeeze(isopl)) 
ncwrite(file_out,'time',anno') 
ncwrite(file_out,'wglvb_occitens_isolev',w_isop) 

% Global Attributes
ncwriteatt(file_out,'/','description',                  'OCCITENS vertical velocities');
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
ncwriteatt(file_out,'wglvb_occitens_isolev',  'axis',         'X,Y,X,T')
ncwriteatt(file_out,'wglvb_occitens_isolev',  'long_name',    'beta-plane geostrophic vertical velocities')
ncwriteatt(file_out,'wglvb_occitens_isolev',  'units','m/s')

% %depth isopycnal surfaces
% ncwriteatt(file_out,'h_isolev',  'axis',         'X,Y,X,T')
% ncwriteatt(file_out,'h_isolev',  'long_name',    'depth isopycnal surfaces')
% ncwriteatt(file_out,'h_isolev',  'units','m')

%% Filter variables =====================================================
% Yearly isopycnal surface:
nfilt       = 5*4/2;
anno        = 1993:2015;

for an = 1:length(anno)
    disp(['Year: ',num2str(anno(an))])

    start = [1,1,1,an];
    count = [Inf,Inf,Inf,1];

    % OLIV3 loading:

    filo    = 'E:\OLIV3_figures\OCCITENS\occitens_wglvb_glob_025_annual_isolevm.nc';
    w_an    = ncread(filo,'wglvb_occitens_isolev',start,count);
    
    dimlon = size(w_an,1); dimlat = size(w_an,2); dimdep = size(w_an,3);
    
    for k = 1:dimdep
        wk              = squeeze(w_an(:,:,k));
        wf(:,:,k,an)    = smooth2a(wk,nfilt);
    end
end
%% Create .nc file ----------------------------------------------------
file_out = 'E:\OLIV3_figures\OCCITENS\occitens_wglvb_glob_025_annual_isolevm_filtr.nc';
nccreate(file_out,'longitude',...
         'Dimensions', {'x',dimlon,'y',dimlat},'Format','netcdf4','Datatype','single')
nccreate(file_out,'latitude',...
         'Dimensions', {'x',dimlon,'y',dimlat},'Format','netcdf4','Datatype','single')
nccreate(file_out,'isolev',...
         'Dimensions', {'lev',length(isopl)},'Format','netcdf4','Datatype','single')
nccreate(file_out,'time',...
         'Dimensions', {'time',length(anno)},'Format','netcdf4','Datatype','single')
nccreate(file_out,'wglvb_occitens_isolev_filtr',...
         'Dimensions', {'x',dimlon,'y',dimlat,'lev',length(isopl),'time',length(anno)},'Format','netcdf4','Datatype','double')

% Write variables
ncwrite(file_out,'longitude',X)
ncwrite(file_out,'latitude',Y)    
ncwrite(file_out,'isolev',squeeze(isopl)) 
ncwrite(file_out,'time',anno') 
ncwrite(file_out,'wglvb_occitens_isolev_filtr',wf) 

% Global Attributes
ncwriteatt(file_out,'/','description',                  'OCCITENS vertical velocities');
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
ncwriteatt(file_out,'wglvb_occitens_isolev_filtr',  'axis',         'X,Y,X,T')
ncwriteatt(file_out,'wglvb_occitens_isolev_filtr',  'long_name',    'beta-plane geostrophic vertical velocities 5 deg filtered')
ncwriteatt(file_out,'wglvb_occitens_isolev_filtr',  'units','m/s')
save('E:\\OLIV3_figures\running\wglvb_filtered.mat','-v7.3','wf')

%%
for an = 1:length(anno)
    disp(['Year: ',num2str(anno(an))])

    start = [1,1,1,an];
    count = [Inf,Inf,Inf,1];

    % OLIV3 loading:

    filo    = 'E:\OLIV3_figures\OCCITENS\occitens_glob_025_annual_isolevm.nc';
    w_an    = ncread(filo,'w_occitens_isolev',start,count);
    
    dimlon = size(w_an,1); dimlat = size(w_an,2); dimdep = size(w_an,3);
    
    for k = 1:dimdep
        wk              = squeeze(w_an(:,:,k));
        wft(:,:,k,an)    = smooth2a(wk,nfilt);
    end
end
% save('E:\\OLIV3_figures\running\wtot_filtered.mat',wft)
file_out = 'E:\OLIV3_figures\OCCITENS\occitens_glob_025_annual_isolevm_filtr.nc';
create_ncfile(file_out, X, Y, isopl, anno, wft, 'w_occitens_isolev_filtr', 'occitens vertical velocities 5 deg filtered');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               METRICS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Correlation coefficient
file_tot    = 'E:\OLIV3_figures\OCCITENS\occitens_glob_025_annual_isolevm_filtr.nc';
file_glvb   = 'E:\OLIV3_figures\OCCITENS\occitens_wglvb_glob_025_annual_isolevm_filtr.nc';

for k = 1:length(isopl)
    disp(['Isopycnal level: ',num2str(k)])

    start       = [1,1,k,1];
    count       = [Inf,Inf,1,Inf];

    wtot_filtr  = squeeze(ncread(file_tot,'w_occitens_isolev_filtr',start,count));
    wglvb_filtr = squeeze(ncread(file_glvb,'wglvb_occitens_isolev_filtr',start,count));

    for i = 1:size(X,1)
        for j = 1:size(Y,2)
            [R,P]       = corrcoef(squeeze(wtot_filtr(i,j,:)),squeeze(wglvb_filtr(i,j,:)),'rows','complete','Alpha',0.05);
            rho(i,j,k)  = R(1,2);
        end
    end
end

% Figure
figure
colormap(turbo(20))
contourf(squeeze(rho(:,:,21))',[-1:0.1:1],'Edgecolor','none')
colorbar
caxis([-1 1])

%% Relative error time-mean

for k = 1:length(isopl)
    disp(['Isopycnal level: ',num2str(k)])

    start       = [1,1,k,1];
    count       = [Inf,Inf,1,Inf];

    wtot_filtr  = squeeze(ncread(file_tot,'w_occitens_isolev_filtr',start,count));
    wglvb_filtr = squeeze(ncread(file_glvb,'wglvb_occitens_isolev_filtr',start,count));

    rel_err_wglvb_tot(:,:,k) = 100.*((mean(wtot_filtr,3)-mean(wglvb_filtr,3))./(abs(mean(wtot_filtr,3))));
end

%% Figure
k = 21;
figure
colormap(turbo(11))
contourf(abs(rel_err_wglvb_tot(:,:,k))',[0:10:110],'Edgecolor','none')
colorbar
caxis([0 110])

%% Save variables .mat

save('E:\OLIV3_figures\OCCITENS\rho_w_wglvb_occitens_glob_025_annual_isolevm_filtr.mat','X','Y','isopl','rho')
save('E:\OLIV3_figures\OCCITENS\relerr_w_wglvb_occitens_glob_025_annual_isolevm_filtr.mat','X','Y','isopl','rel_err_wglvb_tot')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   INTERPOLATION REGULAR GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

circshift_lim = 1442-190+1;

% Circshift variables:
lonc        = circshift(X,circshift_lim,1); 
latc        = circshift(Y,circshift_lim,1);

rhoc        = double(circshift(rho,circshift_lim,1));
relerrc     = double(circshift(rel_err_wglvb_tot,circshift_lim,1));

% Reshape coordinate matrix
latc        = double(latc); lonc = double(lonc);
latr        = reshape(latc,[size(latc,1)*size(latc,2) 1]);
lonr        = reshape(lonc,[size(lonc,1)*size(lonc,2) 1]);
  
% 0.25º grid
LONim       = (-180:0.25:180);
LATim       = (-80:0.25:80);

[LONimone,LATimone] = meshgrid(LONim,LATim');
LATimone    = double(LATimone); 
LONimone    = double(LONimone);

LATimoner   = reshape(LATimone,[size(LATimone,1)*size(LATimone,2) 1]);
LONimoner   = reshape(LONimone,[size(LONimone,1)*size(LONimone,2) 1]);

rho_int     = zeros(size(LATimone,1),size(LATimone,2),71);
relerr_int  = zeros(size(LATimone,1),size(LATimone,2),71);

for k = 1:length(isopl)
    disp(['Isopycnal level: ',num2str(k)])
    rhocr           = reshape(rhoc(:,:,k),[size(rhoc,1)*size(rhoc,2) 1]);
    rhocintr        = griddata(latr,lonr,squeeze(rhocr),LATimoner,LONimoner);
    rho_int(:,:,k)  = reshape(rhocintr,[size(LATimone,1) size(LATimone,2)]);

    relerrcr        = reshape(relerrc(:,:,k),[size(relerrc,1)*size(relerrc,2) 1]);
    relerrcintr     = griddata(latr,lonr,squeeze(relerrcr),LATimoner,LONimoner);
    relerr_int(:,:,k)  = reshape(relerrcintr,[size(LATimone,1) size(LATimone,2)]);
end

%%
save('E:\OLIV3_figures\OCCITENS\rho_w_wglvb_occitens_glob_025_annual_isolevm_filtr.mat','LONimone','LATimone','isopl','rho_int')
save('E:\OLIV3_figures\OCCITENS\relerr_w_wglvb_occitens_glob_025_annual_isolevm_filtr.mat','LONimone','LATimone','isopl','relerr_int')

%% Decomposition
anno        = 1993:2015;
filo    = 'C:\Users\yago_\Documents\LOCEAN\Data\GLORYS\glorys_glob_025_annual.nc';
lon    = ncread(filo,'longitude');
lat    = ncread(filo,'latitude');
depth   = ncread(filo,'depth');

% zoom Canary Islands:
lon_min = -30;
lon_max = -5;
lat_min = 5;
lat_max = 40;

lon_sel = lon(lon(:,1)>lon_min & lon(:,1) < lon_max, lat(1,:)> lat_min & lat(1,:) < lat_max);
lat_sel = lat(lon(:,1)>lon_min & lon(:,1) < lon_max, lat(1,:)> lat_min & lat(1,:) < lat_max);

ww = zeros(size(lon_sel,1),size(lat_sel,2),length(depth),length(anno));
for an = 1:length(anno)
    disp(['Year: ',num2str(anno(an))])

    start = [1,1,1,an];
    count = [Inf,Inf,Inf,1];

    % OLIV3 loading:

    filo    = 'C:\Users\yago_\Documents\LOCEAN\Data\GLORYS\glorys_glob_025_annual.nc';
    w_an    = ncread(filo,'w_glorys',start,count);
    
    ww(:,:,:,an) = w_an(lon(:,1)>lon_min & lon(:,1) < lon_max, lat(1,:)> lat_min & lat(1,:) < lat_max,:);

end
% save('E:\\OLIV3_figures\running\wtot_filtered.mat',wft)
file_out = 'C:\Users\yago_\Documents\LOCEAN\Data\GLORYS\occitens_glob_025_annual_canary.nc';
create_ncfile_verlev(file_out, lon_sel, lat_sel, depth, anno, ww, 'w_glorys', 'glorys12v1 (ORCA025 horizontal grid interpolated) vertical velocities');

% zoom Peru
lon_min = -90;
lon_max = -70;
lat_min = -40;
lat_max = -10;

lon_sel = lon(lon(:,1)>lon_min & lon(:,1) < lon_max, lat(1,:)> lat_min & lat(1,:) < lat_max);
lat_sel = lat(lon(:,1)>lon_min & lon(:,1) < lon_max, lat(1,:)> lat_min & lat(1,:) < lat_max);

ww = zeros(size(lon_sel,1),size(lat_sel,2),length(depth),length(anno));
for an = 1:length(anno)
    disp(['Year: ',num2str(anno(an))])

    start = [1,1,1,an];
    count = [Inf,Inf,Inf,1];

    % OLIV3 loading:

    filo    = 'C:\Users\yago_\Documents\LOCEAN\Data\GLORYS\glorys_glob_025_annual.nc';
    w_an    = ncread(filo,'w_glorys',start,count);
    
    ww(:,:,:,an) = w_an(lon(:,1)>lon_min & lon(:,1) < lon_max, lat(1,:)> lat_min & lat(1,:) < lat_max,:);

end
% save('E:\\OLIV3_figures\running\wtot_filtered.mat',wft)
file_out = 'C:\Users\yago_\Documents\LOCEAN\Data\GLORYS\occitens_glob_025_annual_peru.nc';
create_ncfile_verlev(file_out, lon_sel, lat_sel, depth, anno, ww, 'w_glorys', 'glorys12v1 (ORCA025 horizontal grid interpolated) vertical velocities');
