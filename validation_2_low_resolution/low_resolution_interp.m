% 5 x 5 box analysis
clear all; close all; clc;
% Mean
% Correlation 
% Variance
% Sign change depth

% All variables over [-180:180] grid
%               name            hor. cov.       time. cov
% OCCITENS:     wo_tt/wol_tt    [120:120]       [1993:2019]/[1993:2015]
% OMEGA3D:      wom_tt          [-180:180]      [1993:2019]
% OLIV3:        wa_tt           [0:360]         [1993:2019]
% ECCO:         we_tt           [-37:-37]
% GLORYS:       

% --------------------------------------------------------------------
% BOX 5x5 definition
lon_box5    = [-180:5:180]; lat_box5 = [-80:5:80];
anno23      = [1993:2015];
ind         = [1:54]; isop = log(ind)./1.26 + 25;
isopl1      = [21:0.25:26]; isop_uni = [isopl1 isop(5:end)];
isopl       = isop_uni;
% --------------------------------------------------------------------

%% OCCITENS -----------------------------------------------------------
anno = 1993:2015;

disp('OCCITENS')
path_occ = 'E:\\OLIV3_figures\\OCCITENS\\';

% OCCITENS w total
file_occ    = strcat(path_occ,'occitens_glob_025_annual_isolevm.nc');
lon         = ncread(file_occ,'longitude');
lat         = ncread(file_occ,'latitude');
for an = 1 : length(anno)
    disp(['Year: ',num2str(anno(an))])
    start = [1, 1, 1, an];
    count = [Inf, Inf, Inf, 1];
    
    w1           = ncread(file_occ,'w_occitens_isolev',start,count);
    %h_isop      = ncread(file_occ,'h_isolev',start,count);
  
    % Shift longitude
    circshift_lim   = 1442-430;
    lon         = circshift(lon,circshift_lim,1); 
    lat         = circshift(lat,circshift_lim,1);
    w           = circshift(w,circshift_lim,1);
    lon1        = lon(:,300); lat1 = lat(300,:);
    
    wb_occitens(:,:,:,an) = box_interp_nt(lon1,lat1,isopl,w,lon_box5,lat_box5);
end

save('C:\\Users\\yago_\\Documents\\LOCEAN\\Data\\OCCITENS\\occitens_glob_5_annual_isolevm.mat','wb_occitens')

%% Mixed layer mask OCCITENS
h_isop = ncread('C:\\Users\\yago_\\Documents\\LOCEAN\\Data\\OCCITENS\\occitens_glob_h_isolevm.nc','h_isolev'); % Load isopycnal depths
mmld = nan([size(lon) length(anno)]);
ii = 0;
for an = anno
    ii = ii+1;
    filemld = strcat('C:\\Users\\yago_\\Documents\\LOCEAN\Data\OCCITENS\\mld_annual_means\\ORCA025.L75_annual_mld_',num2str(an),'.nc');
    mmld(:,:,ii) = ncread(filemld,'somxl010');
end
mld = max(mmld,[],3); % Maximum of the whole period

mld_cont    = mld_mask(lon,lat,isopl,mld,h_isop);
mld_cont    = circshift(mld_cont,circshift_lim,1); 

mld_cont    = mld_mask_0(lon,lat,isopl,mld,h_isop);
mld_cont    = circshift(mld_cont,circshift_lim,1); 
mld_cont_b  = box_interp_nt(lon1,lat1,isopl,mld_cont,lon_box5,lat_box5);

mld_cont_bm = mld_cont_b;
mld_cont_bm(mld_cont_bm<0.5) = NaN;
mld_cont_bm(mld_cont_bm>0.5) = 1;

% Save mld mask native grid:

save('C:\\Users\\yago_\\Documents\\LOCEAN\\Data\\OCCITENS\\oliv3_glob_5_mldmask.mat','mld_cont_bm')

% Low resolution interpolation:
wb_occitens_womld = box_interp_woMLD(lon1,lat1,isopl,anno23,w,mld_cont,lon_box5,lat_box5);
save('C:\\Users\\yago_\\Documents\\LOCEAN\\Data\\OCCITENS\\occitens_glob_5_annual_isolevm_woMLD.mat','wb_occitens_womld')

% OCCITENS w glvb --------------------------------------------------------------------
for an = 1 : length(anno)
    disp(['Year: ',num2str(anno(an))])
    start = [1, 1, 1, an];
    count = [Inf, Inf, Inf, 1];
    file_ocl    = strcat(path_occ,'occitens_wglvb_glob_025_annual_isolevm.nc');
    w           = ncread(file_ocl,'wglvb_occitens_isolev',start,count);
    w           = circshift(w,circshift_lim,1);
    wb_glvb(:,:,:,an)     = box_interp_nt(lon1,lat1,isopl,w,lon_box5,lat_box5);
end

save('C:\\Users\\yago_\\Documents\\LOCEAN\\Data\\OCCITENS\\occitens_wglvb_glob_5_annual_isolevm.mat','wb_glvb')

% Interpolation:
for an = 1:length(anno)
    disp(['Year: ',num2str(anno(an))])
    start = [1,1,1,an];
    count = [Inf,Inf,Inf,1];
    file_ocl    = strcat(path_occ,'occitens_wglvb_glob_025_annual_isolevm.nc');
    w           = ncread(file_ocl,'wglvb_occitens_isolev',start,count);
    w           = circshift(w,circshift_lim,1);
    wb_glvb_womld(:,:,:,an) = box_interp_woMLD_nt(lon1,lat1,isopl,w,mld_cont,lon_box5,lat_box5);
end
save('C:\\Users\\yago_\\Documents\\LOCEAN\\Data\\OCCITENS\\occitens_wglvb_glob_5_annual_isolevm_woMLD.mat','wb_glvb_womld')

%% OLIV3 --------------------------------------------------------------
disp('OLIV3')
file_oli    = 'C:\\Users\\yago_\\Documents\\LOCEAN\\Data\\OLIV3\\oliv3_glob_025_ekf_annual_isolevm.nc';

start = [1,1,1,1];
count = [Inf,Inf,Inf,length(anno23)];

w           = ncread(file_oli,'w_oliv3_isolev',start,count);
h_isop      = ncread(file_oli,'h_isolev',start,count);
lon         = ncread(file_oli,'longitude');
lat         = ncread(file_oli,'latitude');

circshift_lim   = length(lon)/2;
lon         = circshift(lon,circshift_lim,1); 
lat         = circshift(lat,circshift_lim,1);
w           = circshift(w,circshift_lim,1);

lon(lon>180) = lon(lon>180) - 360;
lon1        = lon(:,300); lat1 = lat(300,:);


%% Interpolation:
wb_oliv3    = box_interp(lon1,lat1,isopl,anno23,w,lon_box5,lat_box5);
save('C:\\Users\\yago_\\Documents\\LOCEAN\\Data\\OLIV3\\oliv3_glob_5_annual_isolevm.mat','wb_oliv3')

%% Mixed layer mask
% Maximum anual max MLD 
patha = 'E:\\ARMOR3D\\';
mmlda = nan([size(lon,1) size(lat,2) length(anno)]);
ii = 0;
for an = anno
    disp(['Year: ',num2str(an)])
    ii = ii+1;
    filemlda = strcat(patha,'dataset-armor-3d-rep-yearly_',num2str(an),'.mat');
    load(filemlda)
    mmlda(:,:,ii) = mlotstm;
end
mld = max(mmlda,[],3); % Maximum of the whole 56yr period

mld_cont    = mld_mask(lon,lat,isopl,mld,h_isop);

mld_cont    = circshift(mld_cont,circshift_lim,1);

%%
mld_cont    = mld_mask_0(lon,lat,isopl,mld,h_isop);
mld_cont    = circshift(mld_cont,circshift_lim,1); 
mld_cont_b = box_interp_nt(lon1,lat1,isopl,mld_cont,lon_box5,lat_box5);
%%
mld_cont_bm = mld_cont_b;
mld_cont_bm(mld_cont_bm<0.5) = NaN;
mld_cont_bm(mld_cont_bm>0.5) = 1;
save('C:\\Users\\yago_\\Documents\\LOCEAN\\Data\\OCCITENS\\occitens_glob_5_mldmask.mat','mld_cont_bm')
%% Interpolation:
wb_oliv3_womld = box_interp_woMLD(lon1,lat1,isopl,anno23,w,mld_cont,lon_box5,lat_box5);
save('C:\\Users\\yago_\\Documents\\LOCEAN\\Data\\OLIV3\\oliv3_glob_5_annual_isolevm_woMLD.mat','wb_oliv3_womld')

%% OMEGA3D
disp('OMEGA3D')
file_ome    = 'C:\\Users\\yago_\\Documents\\LOCEAN\\Data\\OMEGA3D\\omega3d_glob_025_annual_isolevm.nc';

start = [1,1,1,1];
count = [Inf,Inf,Inf,length(anno23)];

w           = ncread(file_ome,'w_omega3d_isolev',start,count);
lon         = ncread(file_ome,'longitude');
lat         = ncread(file_ome,'latitude');
lon1        = lon(:,300); lat1 = lat(300,:);
%%
wb_omega3d  = box_interp(lon1,lat1,isopl,anno23,w,lon_box5,lat_box5);
save('C:\\Users\\yago_\\Documents\\LOCEAN\\Data\\OMEGA3D\\omega3d_glob_5_annual_isolevm.mat','wb_omega3d')
%%
% Mixed layer mask
% Maximum anual max MLD 
patha = 'E:\\ARMOR3D\\';
mmlda = nan([length(lon) length(lat) length(anno)]);
ii = 0;
for an = anno
    disp(['Year: ',num2str(an)])
    ii = ii+1;
    filemlda = strcat(patha,'dataset-armor-3d-rep-yearly_',num2str(an),'.mat');
    load(filemlda)
    mmlda(:,:,ii) = mlotstm;
end
mld = max(mmlda,[],3); % Maximum of the whole 56yr period

mld_cont    = mld_mask(lon,lat,isopl,mld,h_isop);

%% Interpolation:
wb_omega3d_womld = box_interp_woMLD(lon1,lat1,isopl,anno23,w,mld_cont,lon_box5,lat_box5);
save('C:\\Users\\yago_\\Documents\\LOCEAN\\Data\\OMEGA3D\\omega3d_glob_5_annual_isolevm_woMLD.mat','wb_omega3d_womld')

%% ECCO
disp('ECCO')
file_ecc    = 'C:\\Users\\yago_\\Documents\\LOCEAN\\Data\\ECCO\\ecco_glob_1_annual_isolevm.nc';

start = [1,1,1,1];
count = [Inf,Inf,Inf,length(anno23)];

w           = ncread(file_ecc,'w_ecco_isolev',start,count);
h_isop      = ncread(file_ecc,'h_isolev',start,count);
lon         = ncread(file_ecc,'longitude');
lat         = ncread(file_ecc,'latitude');
%%

wb_ecco     = box_interp(lon,lat,isopl,anno23,w,lon_box5,lat_box5);
save('C:\\Users\\yago_\\Documents\\LOCEAN\\Data\\ECCO\\ecco_glob_5_annual_isolevm.mat','wb_ecco')

%% Mixed layer mask
% Maximum annual max MLD 
mmld = nan([size(lon) length(anno)]);
ii = 0;
for an = anno
    ii = ii+1;
    mmmld = nan([size(lon) 12]);
    for m = 1:12

        filemld = sprintf('E:\\ECCO\\MXLDEPTH\\%d\\MXLDEPTH_%d_%02d.nc',an,an,m);
        mmmld(:,:,m) = ecco_reshape2d(ncread(filemld,'MXLDEPTH'));
    end
    mmld(:,:,ii) = max(mmmld,[],3);
end
mld = max(mmld,[],3); % Maximum of the whole period
%%
mld_cont    = mld_mask(lon,lat,isopl,mld,abs(h_isop));
%%
mld_cont    = mld_mask_0(lon,lat,isopl,mld,abs(h_isop));
mld_cont_b = box_interp_nt(lon,lat,isopl,mld_cont,lon_box5,lat_box5);
%%
mld_cont_bm = mld_cont_b;
mld_cont_bm(mld_cont_bm<0.5) = NaN;
mld_cont_bm(mld_cont_bm>0.5) = 1;
save('C:\\Users\\yago_\\Documents\\LOCEAN\\Data\\OCCITENS\\ecco_glob_5_mldmask.mat','mld_cont_bm')
%%
wb_ecco_womld = box_interp_woMLD(lon,lat,isopl,anno23,w,mld_cont,lon_box5,lat_box5);
save('C:\\Users\\yago_\\Documents\\LOCEAN\\Data\\ECCO\\ecco_glob_5_annual_isolevm_woMLD.mat','wb_ecco_womld')

%% GLORYS
disp('GLORYS')
file_glo    = 'C:\\Users\\yago_\\Documents\\LOCEAN\\Data\\GLORYS\\glorys_glob_025_annual_isolevm.nc';

start = [1,1,1,1];
count = [Inf,Inf,Inf,length(anno23)];

w           = ncread(file_glo,'w_glorys_isolev',start,count);
lon         = ncread(file_glo,'longitude');
lat         = ncread(file_glo,'latitude');

for tt = 1:length(anno23)
    for kk = 1:size(w,3)
        wkt = w(:,:,kk,tt);
        wkt(abs(wkt)>1) = NaN;
        w(:,:,kk,tt) = wkt;
    end
end

% Shift longitude
circshift_lim   = 1442-430;
lon         = circshift(lon,circshift_lim,1); 
lat         = circshift(lat,circshift_lim,1);
w           = circshift(w,circshift_lim,1);
lon1        = lon(:,300); lat1 = lat(300,:);
%%
wb_glorys   = box_interp(lon1,lat1,isopl,anno23,w,lon_box5,lat_box5);
save('C:\\Users\\yago_\\Documents\\LOCEAN\\Data\\GLORYS\\glorys_glob_5_annual_isolevm.mat','wb_glorys')

%% Mixed layer mask
% Maximum annual max MLD 
h_isop = ncread('C:\\Users\\yago_\\Documents\\LOCEAN\\Data\\OCCITENS\\occitens_glob_h_isolevm.nc','h_isolev');
mmld        = nan([size(lon) length(anno23)]);
ii = 0;
for an = anno23
    ii = ii+1;
    filemld = strcat('C:\\Users\\yago_\\Documents\\LOCEAN\Data\OCCITENS\\mld_annual_means\\ORCA025.L75_annual_mld_',num2str(an),'.nc');
    mmld(:,:,ii) = ncread(filemld,'somxl010');
end
mld = max(mmld,[],3); % Maximum of the whole period

mld_cont    = mld_mask(lon,lat,isopl,mld,h_isop);
mld_cont    = circshift(mld_cont,circshift_lim,1); 

% Interpolation:
wb_glorys_womld = box_interp_woMLD(lon1,lat1,isopl,anno23,w,mld_cont,lon_box5,lat_box5);
save('C:\\Users\\yago_\\Documents\\LOCEAN\\Data\\GLORYS\\glorys_glob_5_annual_isolevm_woMLD.mat','wb_glorys_womld')

%% Merging all datasets:

% Including MLD:
wb(:,:,:,:,1)       = cell2mat(struct2cell(load('C:\\Users\\yago_\\Documents\\LOCEAN\\Data\\OLIV3\\oliv3_glob_5_annual_isolevm.mat')));
wb(:,:,:,:,2)       = cell2mat(struct2cell(load('C:\\Users\\yago_\\Documents\\LOCEAN\\Data\\OMEGA3D\\omega3d_glob_5_annual_isolevm.mat')));
wb(:,:,:,:,3)       = cell2mat(struct2cell(load('C:\\Users\\yago_\\Documents\\LOCEAN\\Data\\OCCITENS\\occitens_glob_5_annual_isolevm.mat')));
wb(:,:,:,:,4)       = cell2mat(struct2cell(load('C:\\Users\\yago_\\Documents\\LOCEAN\\Data\\GLORYS\\glorys_glob_5_annual_isolevm.mat')));
wb(:,:,:,:,5)       = cell2mat(struct2cell(load('C:\\Users\\yago_\\Documents\\LOCEAN\\Data\\OCCITENS\\occitens_wglvb_glob_5_annual_isolevm.mat')));
wb(:,:,:,:,6)       = cell2mat(struct2cell(load('C:\\Users\\yago_\\Documents\\LOCEAN\\Data\\ECCO\\ecco_glob_5_annual_isolevm.mat')));
wb(:,:,:,1:end-1,6) = wb(:,:,:,2:end,6);
wb(:,:,:,end,6)     = NaN;
%% Excluding MLD:
wb_wom(:,:,:,:,1)   = cell2mat(struct2cell(load('C:\\Users\\yago_\\Documents\\LOCEAN\\Data\\OLIV3\\oliv3_glob_5_annual_isolevm_woMLD.mat')));
wb_wom(:,:,:,:,2)   = cell2mat(struct2cell(load('C:\\Users\\yago_\\Documents\\LOCEAN\\Data\\OMEGA3D\\omega3d_glob_5_annual_isolevm_woMLD.mat')));
wb_wom(:,:,:,:,3)   = cell2mat(struct2cell(load('C:\\Users\\yago_\\Documents\\LOCEAN\\Data\\OCCITENS\\occitens_glob_5_annual_isolevm_woMLD.mat')));
wb_wom(:,:,:,:,4)   = cell2mat(struct2cell(load('C:\\Users\\yago_\\Documents\\LOCEAN\\Data\\GLORYS\\glorys_glob_5_annual_isolevm_woMLD.mat')));
wb_wom(:,:,:,:,5)   = cell2mat(struct2cell(load('C:\\Users\\yago_\\Documents\\LOCEAN\\Data\\OCCITENS\\occitens_wglvb_glob_5_annual_isolevm_woMLD.mat')));
wb_wom(:,:,:,:,6)   = cell2mat(struct2cell(load('C:\\Users\\yago_\\Documents\\LOCEAN\\Data\\ECCO\\ecco_glob_5_annual_isolevm_woMLD.mat')));
