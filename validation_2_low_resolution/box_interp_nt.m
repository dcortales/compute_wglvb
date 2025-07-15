function w_box = box_interp_nt(lon,lat,depth,w,lon_box,lat_box)

    w_box = nan([length(lon_box) length(lat_box) length(depth)]);

    for k = 1:length(depth)
        disp(['- Isopycnal: ',num2str(depth(k))])
        wtk = squeeze(w(:,:,k));

        for i = 1:length(lon_box)-1
            lon_mask   = nan(size(lon));
            lon_mask(lon>lon_box(i) & lon< lon_box(i+1)) = 1;
            wtki       = wtk.*lon_mask;

            for j = 1:length(lat_box)-1

                lat_mask    = nan(size(lat));
                lat_mask(lat>lat_box(j) & lat< lat_box(j+1)) = 1;
                wtkij       = wtki.*lat_mask;
                w_box(i,j,k)   = mean(wtkij(:),'omitmissing');
            end
        end
    end
end