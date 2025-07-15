function mm = mld_mask_0(lon,lat,isopl,mld,h_isop) 
    
    % Dimensions:
    dimlon = size(lon,1); 
    dimlat = size(lat,2);
    dimdep = length(isopl);

    % Initialise mask:
    mm = zeros([dimlon,dimlat,dimdep]);

    % Mask for each grid point:
    for i = 1:dimlon
        for j = 1:dimlat
            for k = 1:dimdep
                if h_isop(i,j,k) < mld(i,j)
                    mm(i,j,k) = 1;
                elseif h_isop(i,j,k) > mld(i,j)
                    mm(i,j,k) = 0;
                elseif isnan(h_isop(i,j,k)) == 1
                    mm(i,j,k) = NaN;
                end
            end
        end
    end
end