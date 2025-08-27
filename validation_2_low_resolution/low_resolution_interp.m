% 5 x 5 box analysis
clear all; close all; clc;

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
anno        = [1993:2015];
ind         = [1:54]; isop = log(ind)./1.26 + 25;
isopl1      = [21:0.25:26]; isop_uni = [isopl1 isop(5:end)];
isopl       = isop_uni;
% --------------------------------------------------------------------

%% OCCITENS -----------------------------------------------------------

disp('OCCITENS')
path_occ    = '...\OCCITENS\';

% OCCITENS w total
file_occ    = strcat(path_occ,'occitens_glob_025_annual_isolevm.nc');
lon         = ncread(file_occ,'longitude');
lat         = ncread(file_occ,'latitude');
for an = 1 : length(anno)
    disp(['Year: ',num2str(anno(an))])
    start = [1, 1, 1, an];
    count = [Inf, Inf, Inf, 1];
    
    w1           = ncread(file_occ,'w_occitens_isolev',start,count);
  
    % Shift longitude
    circshift_lim   = 1442-430;
    lon         = circshift(lon,circshift_lim,1); 
    lat         = circshift(lat,circshift_lim,1);
    w           = circshift(w,circshift_lim,1);
    lon1        = lon(:,300); lat1 = lat(300,:);
    
    wb_occitens(:,:,:,an) = box_interp_nt(lon1,lat1,isopl,w,lon_box5,lat_box5);
end

save('...\OCCITENS\occitens_glob_5_annual_isolevm.mat','wb_occitens')

%% Mixed layer mask OCCITENS
h_isop = ncread('...\OCCITENS\occitens_glob_h_isolevm.nc','h_isolev'); % Load isopycnal depths computed in validation_1_annual_means
mmld = nan([size(lon) length(anno)]);
ii = 0;
for an = anno
    ii = ii+1;
    filemld = strcat('...\OCCITENS\\mld_annual_means\ORCA025.L75_annual_mld_',num2str(an),'.nc'); % Load MLD data
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

save('...\OCCITENS\oliv3_glob_5_mldmask.mat','mld_cont_bm')

% Low resolution interpolation:

wb_occitens_womld = box_interp_woMLD(lon1,lat1,isopl,anno23,w,mld_cont,lon_box5,lat_box5);
save('...\OCCITENS\occitens_glob_5_annual_isolevm_woMLD.mat','wb_occitens_womld')

%% OCCITENS w glvb --------------------------------------------------------------------
for an = 1 : length(anno)
    disp(['Year: ',num2str(anno(an))])
    start = [1, 1, 1, an];
    count = [Inf, Inf, Inf, 1];
    file_ocl    = strcat(path_occ,'occitens_wglvb_glob_025_annual_isolevm.nc');
    w           = ncread(file_ocl,'wglvb_occitens_isolev',start,count);
    w           = circshift(w,circshift_lim,1);
    wb_glvb(:,:,:,an)     = box_interp_nt(lon1,lat1,isopl,w,lon_box5,lat_box5);
end

save('...\OCCITENS\occitens_wglvb_glob_5_annual_isolevm.mat','wb_glvb')

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
save('...\OCCITENS\occitens_wglvb_glob_5_annual_isolevm_woMLD.mat','wb_glvb_womld')

%% OLIV3 --------------------------------------------------------------
file_oli    = '...\OLIV3\oliv3_glob_025_ekf_annual_isolevm.nc';

w           = ncread(file_oli,'w_oliv3_isolev');     # Load OLIV3 velocities
h_isop      = ncread(file_oli,'h_isolev');           # Load depth isopycnal levels
lon         = ncread(file_oli,'longitude');
lat         = ncread(file_oli,'latitude');

circshift_lim   = length(lon)/2;
lon         = circshift(lon,circshift_lim,1); 
lat         = circshift(lat,circshift_lim,1);
w           = circshift(w,circshift_lim,1);

lon(lon>180) = lon(lon>180) - 360;
lon1        = lon(:,300); lat1 = lat(300,:);

% Low resolution interpolation:
wb_oliv3    = box_interp(lon1,lat1,isopl,anno23,w,lon_box5,lat_box5);
save('...\OLIV3\oliv3_glob_5_annual_isolevm.mat','wb_oliv3')

% Mixed layer mask:
patha = '...\ARMOR3D\';
mmlda = nan([size(lon,1) size(lat,2) length(anno)]);
ii = 0;
for an = anno
    disp(['Year: ',num2str(an)])
    ii = ii+1;
    filemlda = strcat(patha,'dataset-armor-3d-rep-yearly_',num2str(an),'.mat');
    load(filemlda)
    mmlda(:,:,ii) = mlotstm;
end
mld         = max(mmlda,[],3); % Maximum of the whole 56yr period
mld_cont    = mld_mask_0(lon,lat,isopl,mld,h_isop);
mld_cont    = circshift(mld_cont,circshift_lim,1); 
mld_cont_b  = box_interp_nt(lon1,lat1,isopl,mld_cont,lon_box5,lat_box5);

mld_cont_bm = mld_cont_b;
mld_cont_bm(mld_cont_bm<0.5) = NaN;
mld_cont_bm(mld_cont_bm>0.5) = 1;
save('...\ARMOR3D\\occitens_glob_5_mldmask.mat','mld_cont_bm')

% Low resolution interpolation without Mixed-Layer:
wb_oliv3_womld = box_interp_woMLD(lon1,lat1,isopl,anno23,w,mld_cont,lon_box5,lat_box5);
save('...\OLIV3\oliv3_glob_5_annual_isolevm_woMLD.mat','wb_oliv3_womld')

%% OMEGA3D -------------------------------------------------------------------
disp('OMEGA3D')
file_ome    = '...\OMEGA3D\omega3d_glob_025_annual_isolevm.nc';

start = [1,1,1,1];
count = [Inf,Inf,Inf,length(anno23)];

w           = ncread(file_ome,'w_omega3d_isolev',start,count);
lon         = ncread(file_ome,'longitude');
lat         = ncread(file_ome,'latitude');
lon1        = lon(:,300); lat1 = lat(300,:);
%%
wb_omega3d  = box_interp(lon1,lat1,isopl,anno23,w,lon_box5,lat_box5);
save('...\OMEGA3D\omega3d_glob_5_annual_isolevm.mat','wb_omega3d')
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
save('...\OMEGA3D\omega3d_glob_5_annual_isolevm_woMLD.mat','wb_omega3d_womld')

%% ECCO -----------------------------------------------------------
disp('ECCO')
file_ecc    = '...\ECCO\ecco_glob_1_annual_isolevm.nc';

start = [1,1,1,1];
count = [Inf,Inf,Inf,length(anno23)];

w           = ncread(file_ecc,'w_ecco_isolev',start,count);
h_isop      = ncread(file_ecc,'h_isolev',start,count);
lon         = ncread(file_ecc,'longitude');
lat         = ncread(file_ecc,'latitude');
%%

wb_ecco     = box_interp(lon,lat,isopl,anno23,w,lon_box5,lat_box5);
save('...\ECCO\ecco_glob_5_annual_isolevm.mat','wb_ecco')

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
save('...\ECCO\ecco_glob_5_mldmask.mat','mld_cont_bm')
%%
wb_ecco_womld = box_interp_woMLD(lon,lat,isopl,anno23,w,mld_cont,lon_box5,lat_box5);
save('...\ECCO\ecco_glob_5_annual_isolevm_woMLD.mat','wb_ecco_womld')

%% GLORYS ---------------------------------------------------------------
disp('GLORYS')
file_glo    = '...\GLORYS\glorys_glob_025_annual_isolevm.nc';

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
save('...\GLORYS\glorys_glob_5_annual_isolevm.mat','wb_glorys')

%% Mixed layer mask
% Maximum annual max MLD 
h_isop = ncread('...\OCCITENS\occitens_glob_h_isolevm.nc','h_isolev');
mmld        = nan([size(lon) length(anno23)]);
ii = 0;
for an = anno23
    ii = ii+1;
    filemld = strcat('...\OCCITENS\mld_annual_means\ORCA025.L75_annual_mld_',num2str(an),'.nc');
    mmld(:,:,ii) = ncread(filemld,'somxl010');
end
mld = max(mmld,[],3); % Maximum of the whole period

mld_cont    = mld_mask(lon,lat,isopl,mld,h_isop);
mld_cont    = circshift(mld_cont,circshift_lim,1); 

% Interpolation:
wb_glorys_womld = box_interp_woMLD(lon1,lat1,isopl,anno23,w,mld_cont,lon_box5,lat_box5);
save('...\GLORYS\glorys_glob_5_annual_isolevm_woMLD.mat','wb_glorys_womld')

%% Merging all datasets:

% Including MLD:
wb(:,:,:,:,1)       = cell2mat(struct2cell(load('...\OLIV3\oliv3_glob_5_annual_isolevm.mat')));
wb(:,:,:,:,2)       = cell2mat(struct2cell(load('...\OMEGA3D\omega3d_glob_5_annual_isolevm.mat')));
wb(:,:,:,:,3)       = cell2mat(struct2cell(load('...\OCCITENS\occitens_glob_5_annual_isolevm.mat')));
wb(:,:,:,:,4)       = cell2mat(struct2cell(load('...\GLORYS\glorys_glob_5_annual_isolevm.mat')));
wb(:,:,:,:,5)       = cell2mat(struct2cell(load('...\OCCITENS\occitens_wglvb_glob_5_annual_isolevm.mat')));
wb(:,:,:,:,6)       = cell2mat(struct2cell(load('...\ECCO\ecco_glob_5_annual_isolevm.mat')));
wb(:,:,:,1:end-1,6) = wb(:,:,:,2:end,6);
wb(:,:,:,end,6)     = NaN;
%% Excluding MLD:
wb_wom(:,:,:,:,1)   = cell2mat(struct2cell(load('...\OLIV3\oliv3_glob_5_annual_isolevm_woMLD.mat')));
wb_wom(:,:,:,:,2)   = cell2mat(struct2cell(load('...\OMEGA3D\omega3d_glob_5_annual_isolevm_woMLD.mat')));
wb_wom(:,:,:,:,3)   = cell2mat(struct2cell(load('...\OCCITENS\occitens_glob_5_annual_isolevm_woMLD.mat')));
wb_wom(:,:,:,:,4)   = cell2mat(struct2cell(load('...\GLORYS\glorys_glob_5_annual_isolevm_woMLD.mat')));
wb_wom(:,:,:,:,5)   = cell2mat(struct2cell(load('...\OCCITENS\occitens_wglvb_glob_5_annual_isolevm_woMLD.mat')));
wb_wom(:,:,:,:,6)   = cell2mat(struct2cell(load('...\ECCO\ecco_glob_5_annual_isolevm_woMLD.mat')));


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   INTERCOMPARISON RESULTS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Mean
% Variance
% Correlation coefficient
% RMSE
% vertical structure (vertical slope)

%% Mean -------------------------------------------------------------
wm                      = squeeze(mean(wb,4));
[lat_box5m,lon_box5m]   = meshgrid(lat_box5,lon_box5);

save('...\intercomparison_metrics\mean_w.mat','wm','lon_box5m','lat_box5m')

%% Mean without MLD -------------------------------------------------------------
wm                      = squeeze(mean(wb_wom,4));
[lat_box5m,lon_box5m]   = meshgrid(lat_box5,lon_box5);

save('...\intercomparison_metrics\mean_w_womld.mat','wm','lon_box5m','lat_box5m')

%% Variance ---------------------------------------------------------
wv                      = log10(var(wb,0,4,'omitmissing'));
[lat_box5m,lon_box5m]   = meshgrid(lat_box5,lon_box5);

save('...\intercomparison_metrics\variance_w.mat','wv','lon_box5m','lat_box5m')

%% Correlation coefficient ------------------------------------------

rho     = nan([length(lon_box5) length(lat_box5) 71 6 6]);
pval    = nan([length(lon_box5) length(lat_box5) 71 6 6]);

for k = 1:length(isopl)
    for i = 1:length(lon_box5)
        for j = 1:length(lat_box5)
            for ii = 1:5
                for jj = ii+1:6
                    [R,P]               = corrcoef(squeeze(wb(i,j,k,:,ii)),squeeze(wb(i,j,k,:,jj)),'rows','complete','Alpha',0.05);
                    rho(i,j,k,ii,jj)    = R(1,2);
                    pval(i,j,k,ii,jj)   = P(1,2);
                end
            end
        end
    end
end

% Significance ---------------------------------------------------

rho(rho == 0)   = NaN;
pval(pval == 0) = NaN;
%
alfa            = 0.05;
matsigcorr      = nan(size(pval)); 

for i = 1:length(lon_box5)
    for j = 1:length(lat_box5)
        for ii = 1:5
            for jj = 1:6
                for k = 1:length(isopl)
                    if pval(i,j,k,ii,jj) <= alfa
                        matsigcorr(i,j,k,ii,jj) = 1;
                    else
                        matsigcorr(i,j,k,ii,jj) = NaN;
                    end 
                end
            end
        end
    end
end

% RMSE ------------------------------------------------------------

RMSE = nan([length(lon_box5) length(lat_box5) 71 6 6]);
 
for ii = 1:5
    for jj = ii+1:6
        RMSE(:,:,:,jj,ii) = sqrt(mean((wb(:,:,:,:,jj) - wb(:,:,:,:,ii)).^2,4,'omitmissing'));
    end
end

save('...\intercomparison_metrics\R_sign_RMSE_w.mat','rho','matsigcorr','RMSE','lon_box5m','lat_box5m')

%% Vertical gradient ------------------------------------------------- 
% Slope fixed (25-27)
lower_ind = 22; %27 sigma

upper_ind       = ones([length(lon_box5) length(lat_box5) 6]).*29;
m_55            = nan([length(lon_box5) length(lat_box5) 6]); 

for i = 1:length(lon_box5)
    for j = 1:length(lat_box5)
        for ij = 1:6
            ao = 0;
            for k = 1:length(isopl)-1
                if isnan(wm(i,j,k,ij)) == 0 && ao == 0
                    ao = 1;
                    upper_ind(i,j,ij) = k+1;
                end
            end
            if ao == 1
                m_55(i,j,ij) = (abs(wm(i,j,lower_ind,ij))-abs(wm(i,j,upper_ind(i,j,ij),ij)))/(isopl(lower_ind)-isopl(upper_ind(i,j,ij)));
            else
                m_55(i,j,ij) = NaN;
            end
        end
    end
end

m_55(abs(m_55)>4*10^(-6))     = NaN;

save('...\intercomparison_metrics\slope_w_55new_26_5.mat','m_55')

