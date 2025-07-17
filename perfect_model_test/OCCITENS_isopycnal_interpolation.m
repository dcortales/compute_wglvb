%% --------------------------------------------------------------------
% OCCITENS w total
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

    filo    = sprintf('C:\\Users\\yago_\\Documents\\LOCEAN\\Data\\OCCITENS\\w_annual_%04d.nc',anno(an));
    w_an    = permute(ncread(filo,'w'),[3 2 1]);

    % Neutral density loading:
    filesn  = sprintf('C:\\Users\\yago_\\Documents\\LOCEAN\\Data\\OCCITENS\\neutral_density\\ORCA025.L75-OCCITENS.003-signtr_%04d.nc',anno(an));
    sigma   = ncread(filesn,'signtr');

    [w_isop(:,:,:,an), h_isop(:,:,:,an)] = isop_interp(w_an,dimlon,dimlat,dimdep,isopl,Z,sigma);
end

% Create .nc file ----------------------------------------------------
file_out = 'C:\\Users\\yago_\\Documents\\LOCEAN\\Data\\OCCITENS\\occitens_glob_025_annual_isolevy.nc';
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

    [w_isop(:,:,:,an), h_isop] = isop_interp(w_an,dimlon,dimlat,dimdep,isopl,Z,sigma);
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
         'Dimensions', {'x',dimlon,'y',dimlat,'lev',length(isopl)},'Format','netcdf4','Datatype','double')

% Write variables
ncwrite(file_out,'longitude',X)
ncwrite(file_out,'latitude',Y)    
ncwrite(file_out,'isolev',squeeze(isopl)) 
ncwrite(file_out,'time',anno') 
ncwrite(file_out,'w_occitens_isolev',w_isop) 
ncwrite(file_out,'h_isolev',w_isop) 

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

%% Create .nc file ----------------------------------------------------
file_out = 'C:\\Users\\yago_\\Documents\\LOCEAN\\Data\\OCCITENS\\occitens_glob_h_isolevm.nc';
nccreate(file_out,'longitude',...
         'Dimensions', {'x',dimlon,'y',dimlat},'Format','netcdf4','Datatype','single')
nccreate(file_out,'latitude',...
         'Dimensions', {'x',dimlon,'y',dimlat},'Format','netcdf4','Datatype','single')
nccreate(file_out,'isolev',...
         'Dimensions', {'lev',length(isopl)},'Format','netcdf4','Datatype','single')
nccreate(file_out,'h_isolev',...
         'Dimensions', {'x',dimlon,'y',dimlat,'lev',length(isopl)},'Format','netcdf4','Datatype','double')

% Write variables
ncwrite(file_out,'longitude',X)
ncwrite(file_out,'latitude',Y)    
ncwrite(file_out,'isolev',squeeze(isopl)) 
ncwrite(file_out,'h_isolev',h_isop) 

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