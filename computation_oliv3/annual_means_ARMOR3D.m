clear all; close all; clc

%% Cargamos los datos de ARMOR3D.
%Hacemos la media para cada anno y luego para todo el periodo

year        = 1993:2019;
file_path   ='E:\ARMOR3D\1993\dataset-armor-3d-rep-weekly_19930106T1200Z_P20201023T1646Z.nc';
X           = ncread(file_path,'longitude'); 
Y           = ncread(file_path,'latitude');
Z           = ncread(file_path,'depth');
dimlon      = length(X); dimlat = length(Y); dimdep = length(Z);

for y = 20:length(year)
    disp(year(y))
    file_path = cat(2,'E:\ARMOR3D\',num2str(year(y)),'\','dataset-armor-3d-rep-weekly_*.nc');

    %ncdisp(file_path)

    filename = dir(file_path);
    numfiles = length(filename);

    for n = 1:numfiles
       filenames_cell(n,1) = {filename(n).name};
    end

    filenames_cell = char(filenames_cell);
    filenames_cell = string(filenames_cell);

    to = 0; so = 0; ugo = 0; vgo = 0;
    clear tom som mlotstm ugom vgom

    for n = 1:numfiles
        %if n ~= 3
        file = strcat('E:\ARMOR3D\',num2str(year(y)),'\',filenames_cell(n));
        disp(file)
        % ncdisp(file)

        to     = to + ncread(file,'to');
        so     = so + ncread(file,'so');
        ugo    = ugo + ncread(file,'ugo');
        vgo    = vgo + ncread(file,'vgo');
        mlotstm(:,:,n) = ncread(file,'mlotst');
        %end
    end

    tom     = to/numfiles;
    som     = so/numfiles;
    ugom    = ugo/numfiles;
    vgom    = vgo/numfiles;
    %mlotstm = mlotst/numfiles;

    f_out = sprintf('C:\\Users\\yago_\\Documents\\LOCEAN\\Data\\ARMOR3D\\annual_means\\dataset-armor-3d-rep-yearly_%04d.nc',year(y));
    %save(filename_data,'tom','som','ugom','vgom','zom','mlotstm');

    %====================================================
    %         Write data into NetCDF files 
    %====================================================
    
    fillValue   = -9999;
    
    % Define variable dimensions
    nccreate(f_out,'lon',...
        'Dimensions',{'x',dimlon},...
        'Format','netcdf4','Datatype','single')
    nccreate(f_out,'lat',...
        'Dimensions',{'y',dimlat},...
        'Format','netcdf4','Datatype','single')
    nccreate(f_out,'depth',...
        'Dimensions',{'z',dimdep},...
        'Format','netcdf4','Datatype','single')
    nccreate(f_out,'time',...
        'Dimensions',{'tt',numfiles},...
        'Format','netcdf4','Datatype','single')
    nccreate(f_out,'t',...
        'Dimensions',{'x',dimlon,'y',dimlat,'z',dimdep},...
        'Format','netcdf4','Datatype','double')
    nccreate(f_out,'s',...
        'Dimensions',{'x',dimlon,'y',dimlat,'z',dimdep},...
        'Format','netcdf4','Datatype','double')
    nccreate(f_out,'ug',...
        'Dimensions',{'x',dimlon,'y',dimlat,'z',dimdep},...
        'Format','netcdf4','Datatype','double')
    nccreate(f_out,'vg',...
        'Dimensions',{'x',dimlon,'y',dimlat,'z',dimdep},...
        'Format','netcdf4','Datatype','double')
    nccreate(f_out,'ml',...
        'Dimensions',{'x',dimlon,'y',dimlat,'tt',numfiles},...
        'Format','netcdf4','Datatype','double')
    
    % Write variables
    ncwrite(f_out,'lon',X)
    ncwrite(f_out,'lat',Y)
    ncwrite(f_out,'depth',Z)
    ncwrite(f_out,'time',1:numfiles)
    ncwrite(f_out,'t',tom)
    ncwrite(f_out,'s',som)
    ncwrite(f_out,'ug',ugom)
    ncwrite(f_out,'vg',vgom)
    ncwrite(f_out,'ml',mlotstm)
    
    % Global attributes
    ncwriteatt(f_out,'/','title','Annual means ARMOR3D data')
    ncwriteatt(f_out,'/','institution','IMEDEA (UIB-CSIC)')
    ncwriteatt(f_out,'/','CreationDate',datestr(now,'yyyy/mm/dd HH:MM:SS'))
    ncwriteatt(f_out,'/','CreatedBy','D. Cortes-Morales (dcortes@imedea.uib-csic.es)')
    
    % Lon
    ncwriteatt(f_out,'lon','long_name','Longitude')
    ncwriteatt(f_out,'lon','units','degrees_east')
    
    % Lat
    ncwriteatt(f_out,'lat','long_name','Latitude')
    ncwriteatt(f_out,'lat','units','degrees_north')
    
    % Temperature
    ncwriteatt(f_out,'t','long_name','Temperature')
    ncwriteatt(f_out,'t','units','degrees_C')

    % Salinity
    ncwriteatt(f_out,'s','long_name','Salinity')
    ncwriteatt(f_out,'s','units','PSU')

    % Ug
    ncwriteatt(f_out,'ug','long_name','Geostrophic zonal velocity')
    ncwriteatt(f_out,'ug','units','m/s')

    % Vg
    ncwriteatt(f_out,'vg','long_name','Geostrophic meridional velocity')
    ncwriteatt(f_out,'vg','units','m/s')

    % Mixed Layer
    ncwriteatt(f_out,'ml','long_name','Mixed Layer depth')
    ncwriteatt(f_out,'ml','units','m')

end

