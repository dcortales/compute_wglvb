%% Perfect Model Test Metrics
% This script computes the metrics for the perfect model test

% Correlation coefficient:

file_tot    = '...\OCCITENS\occitens_glob_025_annual_isolevm_filtr.nc';
file_glvb   = '...\OCCITENS\occitens_wglvb_glob_025_annual_isolevm_filtr.nc';

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

% Relative error time-mean:

for k = 1:length(isopl)
    disp(['Isopycnal level: ',num2str(k)])

    start       = [1,1,k,1];
    count       = [Inf,Inf,1,Inf];

    wtot_filtr  = squeeze(ncread(file_tot,'w_occitens_isolev_filtr',start,count));
    wglvb_filtr = squeeze(ncread(file_glvb,'wglvb_occitens_isolev_filtr',start,count));

    rel_err_wglvb_tot(:,:,k) = 100.*((mean(wtot_filtr,3)-mean(wglvb_filtr,3))./(abs(mean(wtot_filtr,3))));
end

% Save variables .mat:

save('...\OCCITENS\rho_w_wglvb_occitens_glob_025_annual_isolevm_filtr.mat','X','Y','isopl','rho')
save('...\OCCITENS\relerr_w_wglvb_occitens_glob_025_annual_isolevm_filtr.mat','X','Y','isopl','rel_err_wglvb_tot')

% Vertical gradient:

% Load data and average over time

for k = 1:length(isopl)
    disp(['Isopycnal level: ',num2str(k)])

    start       = [1,1,k,1];
    count       = [Inf,Inf,1,Inf];

    wtot_filtr  = squeeze(ncread(file_tot,'w_occitens_isolev_filtr',start,count));
    wglvb_filtr = squeeze(ncread(file_glvb,'wglvb_occitens_isolev_filtr',start,count));
    
    mean_wtot_filtr(:,:,k)  = mean(wtot_filtr,3);
    mean_wglvb_filtr(:,:,k) = mean(wglvb_filtr,3);
end

% MLD mask computation:

h_isop      = ncread('...\OCCITENS\\occitens_glob_h_isolevm.nc','h_isolev');
mmld        = nan([size(X) length(anno)]);

ii = 0;
for an = anno
    ii              = ii+1;
    filemld         = strcat('...\OCCITENS\\mld_annual_means\\ORCA025.L75_annual_mld_',num2str(an),'.nc');
    mmld(:,:,ii)    = ncread(filemld,'somxl010');
end

mld         = max(mmld,[],3); % Maximum of the whole period
mld_cont    = mld_mask(X,Y,isopl,mld,h_isop);

mean_wtot_filtr_masked  = mean_wtot_filtr.*mld_cont;
mean_wglvb_filtr_masked = mean_wglvb_filtr.*mld_cont;

lower_ind = 29; %27 sigma
upper_ind   = ones(size(X)).*29;

% Compute vertical gradient for wtot

m_55_tot    = nan(size(X)); 
wm          = mean_wtot_filtr_masked;

for i = 1:size(X,1)
    for j = 1:size(Y,2)
            ao = 0;
            for k = 1:length(isopl)-1
                if isnan(wm(i,j,k)) == 0 && ao == 0
                    ao = 1;
                    upper_ind(i,j) = k+1;
                end
            end
            if ao == 1
                m_55_tot(i,j) = (abs(wm(i,j,lower_ind))-abs(wm(i,j,upper_ind(i,j))))/(isopl(lower_ind)-isopl(upper_ind(i,j)));
            else
                m_55_tot(i,j) = NaN;
            end
    end
end

% Compute vertical gradient for wglvb

m_55_glvb   = nan(size(X)); 
wm          = mean_wglvb_filtr_masked;

for i = 1:size(X,1)
    for j = 1:size(Y,2)
            ao = 0;
            for k = 1:length(isopl)-1
                if isnan(wm(i,j,k)) == 0 && ao == 0
                    ao = 1;
                    upper_ind(i,j) = k+1;
                end
            end
            if ao == 1
                m_55_glvb(i,j) = (abs(wm(i,j,lower_ind))-abs(wm(i,j,upper_ind(i,j))))/(isopl(lower_ind)-isopl(upper_ind(i,j)));
            else
                m_55_glvb(i,j) = NaN;
            end
    end
end

save('...\intercomparison_metrics\\slope_w_55new.mat','m_55_tot','m_55_glvb','X','Y','isopl')

% Interpolarion to regular 0.25º grid

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

m_55_totcr       = reshape(m_55_totc,[size(m_55_totc,1)*size(m_55_totc,2) 1]);
m_55_totcintr    = griddata(latr,lonr,squeeze(m_55_totcr),LATimoner,LONimoner);
m_55_tot_int     = reshape(m_55_totcintr,[size(LATimone,1) size(LATimone,2)]);

m_55_glvbcr       = reshape(m_55_glvbc,[size(m_55_glvbc,1)*size(m_55_glvbc,2) 1]);
m_55_glvbcintr    = griddata(latr,lonr,squeeze(m_55_glvbcr),LATimoner,LONimoner);
m_55_glvb_int     = reshape(m_55_glvbcintr,[size(LATimone,1) size(LATimone,2)]);

save('...\OCCITENS\rho_w_wglvb_occitens_glob_025_annual_isolevm_filtr.mat','LONimone','LATimone','isopl','rho_int')
save('...\OCCITENS\relerr_w_wglvb_occitens_glob_025_annual_isolevm_filtr.mat','LONimone','LATimone','isopl','relerr_int')
save('...\OCCITENS\vert_grad_w_wglvb_occitens_glob_025_annual_isolevm_filtr.mat','LONimone','LATimone','m_55_tot_int','m_55_glvb_int')

%% Correlation coefficient OCCITENS wek vs wtot

file_path   ='...\OCCITENS\ORCA025.L75-OCCITENS.003_y1993m12.1m_gridW.nc';
X           = ncread(file_path,'nav_lon'); 
Y           = ncread(file_path,'nav_lat');
Z           = ncread(file_path,'depthw');

% Filter Ekman pumping
disp('Ekman pumping filter')
anno        = 1993:2015;
nfilt       = 5*4/2;

for an = 1:length(anno)
    disp(['Year: ',num2str(anno(an))])

    filo    = sprintf('...\OCCITENS\\wek_annual_%04d.nc',anno(an));
    w_an    = ncread(filo,'wek');

    wekf(:,:,an)    = smooth2a(w_an,nfilt);
end

% Isopycnal definition ------------------------------------------------
ind         = [1:54];
isop        = log(ind)./1.26+25;
isopl1      = [21:0.25:26];
isop_uni    = [isopl1 isop(5:end)];
isopl       = isop_uni;
% ----------------------------------------------------------------------

% Correlation coefficient computation
file_tot    = '...\OCCITENS\occitens_glob_025_annual_isolevm_filtr.nc';

for k = 1:length(isopl)
    disp(['Isopycnal level: ',num2str(k)])

    start       = [1,1,k,1];
    count       = [Inf,Inf,1,Inf];

    wtot_filtr  = squeeze(ncread(file_tot,'w_occitens_isolev_filtr',start,count));

    for i = 1:size(X,1)
        for j = 1:size(Y,2)
            [R,P]       = corrcoef(squeeze(wtot_filtr(i,j,:)),squeeze(wekf(j,i,:)),'rows','complete','Alpha',0.05);
            rho(i,j,k)  = R(1,2);
        end
    end
end

circshift_lim = 1442-190+1;

% Circshift variables:
lonc        = circshift(X,circshift_lim,1); 
latc        = circshift(Y,circshift_lim,1);

rhoc        = double(circshift(rho,circshift_lim,1));

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

for k = 19:21
    disp(['Isopycnal level: ',num2str(k)])
    rhocr           = reshape(rhoc(:,:,k),[size(rhoc,1)*size(rhoc,2) 1]);
    rhocintr        = griddata(latr,lonr,squeeze(rhocr),LATimoner,LONimoner);
    rho_int(:,:,k)  = reshape(rhocintr,[size(LATimone,1) size(LATimone,2)]);
end

save('rho_wek_w_occitens_glob_025_annual_isolevm_filtr.mat','LONimone','LATimone','isopl','rho_int')

%% Correlation coefficient OCCITENS wek vs wtot (horizontal levels)

disp('Isopycnal interpolation:')

file_path   ='...\OCCITENS\ORCA025.L75-OCCITENS.003\GRID\1993\ORCA025.L75-OCCITENS.003_y1993m01.1m_gridW.nc';
X           = ncread(file_path,'nav_lon'); 
Y           = ncread(file_path,'nav_lat');
Z           = ncread(file_path,'depthw');

% Isopycnal definition ------------------------------------------------
ind         = [1:54];
isop        = log(ind)./1.26+25;
isopl1      = [21:0.25:26];
isop_uni    = [isopl1 isop(5:end)];
isopl       = isop_uni;
% ----------------------------------------------------------------------

% Filter Ekman pumping
disp('Ekman pumping filter')
anno        = 1993:2015;
nfilt       = 5*4/2;
k = 25;
for an = 1:length(anno)
    disp(['Year: ',num2str(anno(an))])

    start       = [k,1,1];
    count       = [1,Inf,Inf];

    filo    = sprintf('...\OCCITENS\\w_annual_%04d.nc',anno(an));
    w_an    = squeeze(ncread(filo,'w',start,count));

    wtot(:,:,an)    = smooth2a(w_an,nfilt);
end

for an = 1:length(anno)
    disp(['Year: ',num2str(anno(an))])

    filo    = sprintf('...\OCCITENS\\wek_annual_%04d.nc',anno(an));
    w_an    = ncread(filo,'wek');

    wekf(:,:,an)    = smooth2a(w_an,nfilt);
end

for i = 1:size(X,1)
    for j = 1:size(Y,2)
        [R,P]       = corrcoef(squeeze(wtot(j,i,:)),squeeze(wekf(j,i,:)),'rows','complete','Alpha',0.05);
        rho(i,j)  = R(1,2);
    end
end

circshift_lim = 1442-190+1;

% Circshift variables:
lonc        = circshift(X,circshift_lim,1); 
latc        = circshift(Y,circshift_lim,1);

rhoc        = double(circshift(rho,circshift_lim,1));

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

rhocr           = reshape(rhoc,[size(rhoc,1)*size(rhoc,2) 1]);
rhocintr        = griddata(latr,lonr,squeeze(rhocr),LATimoner,LONimoner);
rho_int  = reshape(rhocintr,[size(LATimone,1) size(LATimone,2)]);

save('rho_wek_w_occitens_glob_025_annual_horlev_filtr.mat','LONimone','LATimone','isopl','rho_int')

%% Correlation coefficient OCCITENS wek vs wglvb (horizontal levels)

disp('Isopycnal interpolation:')

file_path   ='...\OCCITENS\ORCA025.L75-OCCITENS.003\GRID\1993\ORCA025.L75-OCCITENS.003_y1993m01.1m_gridW.nc';
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

% Filter Ekman pumping
disp('Ekman pumping filter')
anno        = 1993:2015;
nfilt       = 5*4/2;

k = 25; % 100 m depth
for an = 1:length(anno)
    disp(['Year: ',num2str(anno(an))])

    start       = [k,1,1];
    count       = [1,Inf,Inf];

    filo    = sprintf('...\OCCITENS\\wglvb_annual_%04d.nc',anno(an));
    w_an    = squeeze(ncread(filo,'div',start,count));

    wtot(:,:,an)    = smooth2a(w_an,nfilt);
end

for an = 1:length(anno)
    disp(['Year: ',num2str(anno(an))])

    filo    = sprintf('...\OCCITENS\\wek_annual_%04d.nc',anno(an));
    w_an    = ncread(filo,'wek');

    wekf(:,:,an)    = smooth2a(w_an,nfilt);
end

for i = 1:size(X,1)
    for j = 1:size(Y,2)
        [R,P]       = corrcoef(squeeze(wtot(j,i,:)),squeeze(wekf(j,i,:)),'rows','complete','Alpha',0.05);
        rho(i,j)  = R(1,2);
    end
end

circshift_lim = 1442-190+1;

% Circshift variables:
lonc        = circshift(X,circshift_lim,1); 
latc        = circshift(Y,circshift_lim,1);

rhoc        = double(circshift(rho,circshift_lim,1));

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

rhocr           = reshape(rhoc,[size(rhoc,1)*size(rhoc,2) 1]);
rhocintr        = griddata(latr,lonr,squeeze(rhocr),LATimoner,LONimoner);
rho_int  = reshape(rhocintr,[size(LATimone,1) size(LATimone,2)]);

save('rho_wek_wglvb_occitens_glob_025_annual_horlev_filtr.mat','LONimone','LATimone','isopl','rho_int')
