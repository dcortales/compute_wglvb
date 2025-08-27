% OMEGA3D anual means

clear all; close all; clc

year        = ;
file_path   = '...\OMEGA3D\1993\dataset-omega-3d-rep-weekly_19930106T0000Z_P20200331T0000Z.nc';
X           = ncread(file_path,'lon'); 
Y           = ncread(file_path,'lat');
Z           = ncread(file_path,'depth');
dimlon      = length(X); dimlat = length(Y); dimdep = length(Z);

%% Step 1: Compute annual z-level means

for y = 1:length(year)
    disp(['Year: ', num2str(year(y))])
    file_path   = cat(2,'...\OMEGA3D\',num2str(year(y)),'\','dataset-omega-3d-rep-weekly_*.nc');

    filename    = dir(file_path);
    numfiles    = length(filename);

    for n = 1:numfiles
       filenames_cell(n,1) = {filename(n).name};
    end

    filenames_cell = char(filenames_cell);
    filenames_cell = string(filenames_cell);

    wo = 0;

    for n = 1:numfiles
        file    = strcat('..\OMEGA3D\',num2str(year(y)),'\',filenames_cell(n));
        disp(file)
        wo      = wo + ncread(file,'wo')/86400;
    end
    wom(:,:,:,y) = wo/numfiles;
end

% Write data into NetCDF files:

file_out = '..\OMEGA3D\omega3d_glob_025_annual_verlev.nc';

nccreate(file_out,'longitude',...
         'Dimensions', {'x',dimlon,'y',dimlat},'Format','netcdf4','Datatype','single')
nccreate(file_out,'latitude',...
         'Dimensions', {'x',dimlon,'y',dimlat},'Format','netcdf4','Datatype','single')
nccreate(file_out,'depth',...
         'Dimensions', {'deptht',dimdep},'Format','netcdf4','Datatype','single')
nccreate(file_out,'time',...
         'Dimensions', {'time',length(year)},'Format','netcdf4','Datatype','single')
nccreate(file_out,'w_omega3d',...
         'Dimensions', {'x',dimlon,'y',dimlat,'deptht',dimdep,'time',length(year)},'Format','netcdf4','Datatype','double')

% Write variables
ncwrite(file_out,'longitude',X)
ncwrite(file_out,'latitude',Y)    
ncwrite(file_out,'depth',Z) 
ncwrite(file_out,'time',year') 
ncwrite(file_out,'w_omega3d',wom) 

% Global Attributes
ncwriteatt(file_out,'/','description',          'OMEGA3D vertical velocities');
ncwriteatt(file_out,'/','output_frequency',     '1yr');
    
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
ncwriteatt(file_out,'w_omega3d',  'axis',         'X,Y,X,T')
ncwriteatt(file_out,'w_omega3d',  'long_name',    'quasi-geostrophic OMEGA3D vertical velocities')
ncwriteatt(file_out,'w_omega3d',  'units','m/s')

%% Step 2: Interpolate to yearly isopycnals over ARMOR3D grid

disp('Isopycnal interpolation:')
anno        = 1993:2019;

file_path   ='...\OMEGA3D\1993\dataset-omega-3d-rep-weekly_19930106T0000Z_P20200331T0000Z.nc';
X           = ncread(file_path,'lon'); 
Y           = ncread(file_path,'lat');
Y           = Y(32:end);   % ARMOR3D grid
Z           = ncread(file_path,'depth');
dimlon      = length(X); dimlat = length(Y); dimdep = length(Z);
[LATm,LONm] = meshgrid(double(Y),double(X));

% Isopycnal definition ------------------------------------------------
ind         = [1:54];
isop        = log(ind)./1.26+25;
isopl1      = [21:0.25:26];
isop_uni    = [isopl1 isop(5:end)];
isopl       = isop_uni;
% ----------------------------------------------------------------------

filesn  = sprintf('...\ARMOR3D\neutral_density\dataset-armor-3d-rep-yearly-signtr_%04d.nc',anno(1));
depth   = ncread(filesn,'depth');
lon     = ncread(filesn,'lon');
lonc    = circshift(lon,length(lon)/2);
lonc(1:length(lon)/2) = lonc(1:length(lon)/2)-360;

for an = 1:length(anno)
    disp(['Year: ',num2str(anno(an))])

    start = [1,1,1,an];
    count = [Inf,Inf,Inf,1];

    % OLIV3 loading:

    filo    = '...\OMEGA3D\omega3d_glob_025_annual_verlev.nc';
    w_an    = ncread(filo,'w_omega3d',start,count);

    % Neutral density loading:
    filesn  = sprintf('...\ARMOR3D\neutral_density\dataset-armor-3d-rep-yearly-signtr_%04d.nc',anno(an));
    sigma   = ncread(filesn,'signtr');
    sigmamd = circshift(sigma(:,:,1:41),length(lon)/2,1); % Eliminate ARMOR3D data below 1500m

    sigmamdg = nan(size(w_an(:,32:end,:)));
    for ii = 1:dimlon
        for jj = 1:dimlat
            sigmamdg(ii,jj,:) = interp1(depth(1:41),squeeze(sigmamd(ii,jj,:)),Z);
        end
    end
    
    dimlon = size(w_an,1); dimlat = size(w_an(:,32:end,:),2); dimdep = size(w_an,3);

    [w_isop(:,:,:,an), h_isop(:,:,:,an)] = isop_interp(w_an(:,32:end,:),dimlon,dimlat,dimdep,isopl,Z,sigmamdg);
end

% Write data into NetCDF files:

file_out = '...\OMEGA3D\omega3d_glob_025_annual_isolevy.nc';
nccreate(file_out,'longitude',...
         'Dimensions', {'x',dimlon,'y',dimlat},'Format','netcdf4','Datatype','single')
nccreate(file_out,'latitude',...
         'Dimensions', {'x',dimlon,'y',dimlat},'Format','netcdf4','Datatype','single')
nccreate(file_out,'isolev',...
         'Dimensions', {'lev',length(isopl)},'Format','netcdf4','Datatype','single')
nccreate(file_out,'time',...
         'Dimensions', {'time',length(anno)},'Format','netcdf4','Datatype','single')
nccreate(file_out,'w_omega3d_isolev',...
         'Dimensions', {'x',dimlon,'y',dimlat,'lev',length(isopl),'time',length(anno)},'Format','netcdf4','Datatype','double')
nccreate(file_out,'h_isolev',...
         'Dimensions', {'x',dimlon,'y',dimlat,'lev',length(isopl),'time',length(anno)},'Format','netcdf4','Datatype','double')

% Write variables
ncwrite(file_out,'longitude',LONm)
ncwrite(file_out,'latitude',LATm)    
ncwrite(file_out,'isolev',squeeze(isopl)) 
ncwrite(file_out,'time',anno') 
ncwrite(file_out,'w_omega3d_isolev',w_isop) 

% Global Attributes
ncwriteatt(file_out,'/','description',                  'OMEGA3D vertical velocities');
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
ncwriteatt(file_out,'w_omega3d_isolev',  'axis',         'X,Y,X,T')
ncwriteatt(file_out,'w_omega3d_isolev',  'long_name',    'quasi-geostrophic OMEGA3D vertical velocities')
ncwriteatt(file_out,'w_omega3d_isolev',  'units','m/s')

%depth isopycnal surfaces
ncwriteatt(file_out,'h_isolev',  'axis',         'X,Y,X,T')
ncwriteatt(file_out,'h_isolev',  'long_name',    'ARMOR3D depth isopycnal surfaces')
ncwriteatt(file_out,'h_isolev',  'units','m')
