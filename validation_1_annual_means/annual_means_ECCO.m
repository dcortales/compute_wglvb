%% ECCO
clear all; close all; clc

disp('----------------------------ECCO----------------------------------')
anno        = 1992:2015; month = 1:12;

pathECCO    = '...\\ECCO\\';
pathTT      = string(sprintf('THETA\\%4d\\THETA_%4d_%02d.nc',anno(1),anno(1),month(1)));
TT          = strcat(pathECCO,pathTT);
depthe      = ncread(TT,'Z'); LONe = ncread(TT,'XC'); LATe = ncread(TT,'YC');

lone        = ecco_reshape2d(LONe);
late        = ecco_reshape2d(LATe);

%% Definimos las matrices de las variables.
for an = 1:length(anno)
    w_m = nan(90,90,13,50,12);
    for m = 1:length(month)
        fprintf('%d/%02d\n',anno(an),m)
        file_WVELMASS       = string(sprintf('...\\ECCO\\WVELMASS\\%4d\\WVELMASS_%4d_%02d.nc',anno(an),anno(an),m));
        w_m(:,:,:,:,m)      = squeeze(ncread(file_WVELMASS,'WVELMASS'));

    end
    wme     = squeeze(mean(w_m,5)); 
    we_glo(:,:,:,an)  = ecco_reshape3d(wme);

end

we_glo(we_glo == 0) = NaN;

%% Create .nc file ----------------------------------------------------
file_out = '...\\ECCO\\ecco_glob_1_annual_verlev.nc';
nccreate(file_out,'longitude',...
         'Dimensions', {'x',size(lone,1),'y',size(late,2)},'Format','netcdf4','Datatype','single')
nccreate(file_out,'latitude',...
         'Dimensions', {'x',size(lone,1),'y',size(late,2)},'Format','netcdf4','Datatype','single')
nccreate(file_out,'depth',...
         'Dimensions', {'deptht',size(depthe,1)},'Format','netcdf4','Datatype','single')
nccreate(file_out,'time',...
         'Dimensions', {'time',length(anno)},'Format','netcdf4','Datatype','single')
nccreate(file_out,'w_ecco',...
         'Dimensions', {'x',size(lone,1),'y',size(late,2),'deptht',size(depthe,1),'time',length(anno)},'Format','netcdf4','Datatype','double')

% Write variables
ncwrite(file_out,'longitude',lone)
ncwrite(file_out,'latitude',late)    
ncwrite(file_out,'depth',depthe) 
ncwrite(file_out,'time',anno') 
ncwrite(file_out,'w_ecco',we_glo) 

% Global Attributes
ncwriteatt(file_out,'/','description',          'ECCOv4r4 reanalysis vertical velocities');
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
ncwriteatt(file_out,'w_ecco',  'axis',         'X,Y,X,T')
ncwriteatt(file_out,'w_ecco',  'long_name',    'output vertical velocities')
ncwriteatt(file_out,'w_ecco',  'units','m/s')

%% Temperature-salinity GLOBAL fields

anno        = 2016:2017; month = 1:12;

pathECCO    = '...\\ECCO\\';
pathTT      = string(sprintf('THETA\\%4d\\THETA_%4d_%02d.nc',anno(1),anno(1),month(1)));
TT          = strcat(pathECCO,pathTT);
depthe      = ncread(TT,'Z'); LONe = ncread(TT,'XC'); LATe = ncread(TT,'YC');

lone        = ecco_reshape2d(LONe);
late        = ecco_reshape2d(LATe);

% Definimos las matrices de las variables.
for an = 1:length(anno)
    t_m = nan(90,90,13,50,12);
    s_m = nan(90,90,13,50,12);
    for m = 1:length(month)
        fprintf('%d/%02d\n',anno(an),m)
        file_T          = string(sprintf('...\\ECCO\\THETA\\%4d\\THETA_%4d_%02d.nc',anno(an),anno(an),m));
        t_m  = squeeze(ncread(file_T,'THETA'));
        file_S          = string(sprintf('...\\ECCO\\SALT\\%4d\\SALT_%4d_%02d.nc',anno(an),anno(an),m));
        s_m  = squeeze(ncread(file_S,'SALT'));

    %tme     = squeeze(mean(t_m,5)); 
    %sme     = squeeze(mean(t_m,5)); 
    te_glo(:,:,:)  = ecco_reshape3d(t_m);
    se_glo(:,:,:)  = ecco_reshape3d(s_m);

te_glo(te_glo == 0) = NaN;
se_glo(se_glo == 0) = NaN;

% Create .nc file ----------------------------------------------------
file_out = sprintf('...\\ECCO\\TS\\%4d\\TS_%4d_%02d.nc',anno(an),anno(an),m);
nccreate(file_out,'longitude',...
         'Dimensions', {'x',size(lone,1),'y',size(late,2)},'Format','netcdf4','Datatype','single')
nccreate(file_out,'latitude',...
         'Dimensions', {'x',size(lone,1),'y',size(late,2)},'Format','netcdf4','Datatype','single')
nccreate(file_out,'depth',...
         'Dimensions', {'deptht',size(depthe,1)},'Format','netcdf4','Datatype','single')
nccreate(file_out,'votemper',...
         'Dimensions', {'x',size(lone,1),'y',size(late,2),'deptht',size(depthe,1)},'Format','netcdf4','Datatype','double')
nccreate(file_out,'vosaline',...
         'Dimensions', {'x',size(lone,1),'y',size(late,2),'deptht',size(depthe,1)},'Format','netcdf4','Datatype','double')

% Write variables
ncwrite(file_out,'longitude',lone)
ncwrite(file_out,'latitude',late)    
ncwrite(file_out,'depth',depthe) 
ncwrite(file_out,'votemper',te_glo) 
ncwrite(file_out,'vosaline',se_glo) 

% Global Attributes
ncwriteatt(file_out,'/','description',          'ECCOv4r4 reanalysis temperature and salinity');
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

%temperature and salinity
ncwriteatt(file_out,'votemper',  'axis',         'X,Y,X,T')
ncwriteatt(file_out,'votemper',  'long_name',    'output temperature')
ncwriteatt(file_out,'votemper',  'units','deg Celsius')

ncwriteatt(file_out,'vosaline',  'axis',         'X,Y,X,T')
ncwriteatt(file_out,'vosaline',  'long_name',    'output salinity')
ncwriteatt(file_out,'vosaline',  'units','PSU')
    end
end

%% Isopycnal interpolation

disp('---------------Isopycnal interpolation---------------')
anno        = 1992:2017; month = 1:12;
ind = [1:54];
isop =log(ind)./1.26+25;
isopl1 = [21:0.25:26];
isop_uni = [isopl1 isop(5:end)];
isopl = isop_uni;

% Variables over isopycnals
%signtr = mean(signtr_an,4);
%sigma = signtr;

for an = 1:length(anno)
    disp(['Year: ',num2str(anno(an))])

    start = [1,1,1,an];
    count = [Inf,Inf,Inf,1];

    % ECCO w loading:

    filo    = ...\\ECCO\\ecco_glob_1_annual_verlev.nc';
    w_an    = ncread(filo,'w_ecco',start,count);

    % Neutral density loading:
    filesn  = sprintf('...\\ECCO\\neutral_density\\SIGNTR_GLOB_%04d.nc',anno(an));
    sigma   = ncread(filesn,'signtr');
    
    dimlon = size(w_an,1); dimlat = size(w_an,2); dimdep = size(w_an,3);

    [w_isop(:,:,:,an), h_isop(:,:,:,an)] = isop_interp(w_an,dimlon,dimlat,dimdep,isopl,depthe,sigma);
end

% Create .nc file ----------------------------------------------------
file_out = '...\\ECCO\\ecco_glob_1_annual_isolevy.nc';
nccreate(file_out,'longitude',...
         'Dimensions', {'x',size(lone,1),'y',size(late,2)},'Format','netcdf4','Datatype','single')
nccreate(file_out,'latitude',...
         'Dimensions', {'x',size(lone,1),'y',size(late,2)},'Format','netcdf4','Datatype','single')
nccreate(file_out,'isolev',...
         'Dimensions', {'lev',length(isopl)},'Format','netcdf4','Datatype','single')
nccreate(file_out,'time',...
         'Dimensions', {'time',length(anno)},'Format','netcdf4','Datatype','single')
nccreate(file_out,'w_ecco_isolev',...
         'Dimensions', {'x',size(lone,1),'y',size(late,2),'lev',length(isopl),'time',length(anno)},'Format','netcdf4','Datatype','double')
nccreate(file_out,'h_isolev',...
         'Dimensions', {'x',size(lone,1),'y',size(late,2),'lev',length(isopl),'time',length(anno)},'Format','netcdf4','Datatype','double')

% Write variables
ncwrite(file_out,'longitude',lone)
ncwrite(file_out,'latitude',late)    
ncwrite(file_out,'isolev',squeeze(isopl)) 
ncwrite(file_out,'time',anno') 
ncwrite(file_out,'w_ecco_isolev',w_isop) 
ncwrite(file_out,'h_isolev',h_isop) 

% Global Attributes
ncwriteatt(file_out,'/','description',          'ECCOv4r4 reanalysis vertical velocities');
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
ncwriteatt(file_out,'isolev',    'axis',         'Z')
ncwriteatt(file_out,'isolev',    'long_name',    'Isopycnal levels')
ncwriteatt(file_out,'isolev',    'units',        'kg m^-3')
ncwriteatt(file_out,'isolev',    'positive',     'down')

%time
ncwriteatt(file_out,'time',     'axis',         'T')
ncwriteatt(file_out,'time',     'long_name',    'year')
ncwriteatt(file_out,'time',     'units',        'year')

%beta-plane geostrophic velocities
ncwriteatt(file_out,'w_ecco_isolev',  'axis',         'X,Y,X,T')
ncwriteatt(file_out,'w_ecco_isolev',  'long_name',    'output vertical velocities')
ncwriteatt(file_out,'w_ecco_isolev',  'units','m/s')

ncwriteatt(file_out,'h_isolev',  'axis',         'X,Y,X,T')
ncwriteatt(file_out,'h_isolev',  'long_name',    'depth isopycnal surfaces')
ncwriteatt(file_out,'h_isolev',  'units','m')

%% Isopycnal interpolation

disp('---------------Isopycnal interpolation---------------')
anno        = 1992:2015; month = 1:12;
ind = [1:54];
isop =log(ind)./1.26+25;
isopl1 = [21:0.25:26];
isop_uni = [isopl1 isop(5:end)];
isopl = isop_uni;

% Variables over isopycnals

% Neutral density loading:
filesn  = '...\\ECCO\\neutral_density\\SIGNTR_GLOB_27yr.nc';
sigma   = ncread(filesn,'signtr');

for an = 1%:length(anno)
    disp(['Year: ',num2str(anno(an))])

    start = [1,1,1,an];
    count = [Inf,Inf,Inf,1];

    % ECCO w loading:

    filo    = '...\\ECCO\\ecco_glob_1_annual_verlev.nc';
    w_an    = ncread(filo,'w_ecco',start,count);

    dimlon  = size(w_an,1); dimlat = size(w_an,2); dimdep = size(w_an,3);

    [w_isop(:,:,:,an), h_isop(:,:,:,an)] = isop_interp(w_an,dimlon,dimlat,dimdep,isopl,depthe,sigma);
end

% Create .nc file ----------------------------------------------------
file_out = '...\\ECCO\\ecco_glob_1_annual_isolevm.nc';
nccreate(file_out,'longitude',...
         'Dimensions', {'x',size(lone,1),'y',size(late,2)},'Format','netcdf4','Datatype','single')
nccreate(file_out,'latitude',...
         'Dimensions', {'x',size(lone,1),'y',size(late,2)},'Format','netcdf4','Datatype','single')
nccreate(file_out,'isolev',...
         'Dimensions', {'lev',length(isopl)},'Format','netcdf4','Datatype','single')
nccreate(file_out,'time',...
         'Dimensions', {'time',length(anno)},'Format','netcdf4','Datatype','single')
nccreate(file_out,'w_ecco_isolev',...
         'Dimensions', {'x',size(lone,1),'y',size(late,2),'lev',length(isopl),'time',length(anno)},'Format','netcdf4','Datatype','double')
nccreate(file_out,'h_isolev',...
         'Dimensions', {'x',size(lone,1),'y',size(late,2),'lev',length(isopl),'time',length(anno)},'Format','netcdf4','Datatype','double')

% Write variables
ncwrite(file_out,'longitude',lone)
ncwrite(file_out,'latitude',late)    
ncwrite(file_out,'isolev',squeeze(isopl)) 
ncwrite(file_out,'time',anno') 
ncwrite(file_out,'w_ecco_isolev',w_isop) 
ncwrite(file_out,'h_isolev',h_isop) 

% Global Attributes
ncwriteatt(file_out,'/','description',          'ECCOv4r4 reanalysis vertical velocities');
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
ncwriteatt(file_out,'isolev',    'axis',         'Z')
ncwriteatt(file_out,'isolev',    'long_name',    'Isopycnal levels')
ncwriteatt(file_out,'isolev',    'units',        'kg m^-3')
ncwriteatt(file_out,'isolev',    'positive',     'down')

%time
ncwriteatt(file_out,'time',     'axis',         'T')
ncwriteatt(file_out,'time',     'long_name',    'year')
ncwriteatt(file_out,'time',     'units',        'year')

%beta-plane geostrophic velocities
ncwriteatt(file_out,'w_ecco_isolev',  'axis',         'X,Y,X,T')
ncwriteatt(file_out,'w_ecco_isolev',  'long_name',    'output vertical velocities')
ncwriteatt(file_out,'w_ecco_isolev',  'units','m/s')

ncwriteatt(file_out,'h_isolev',  'axis',         'X,Y,X,T')
ncwriteatt(file_out,'h_isolev',  'long_name',    'depth isopycnal surfaces')
ncwriteatt(file_out,'h_isolev',  'units','m')
