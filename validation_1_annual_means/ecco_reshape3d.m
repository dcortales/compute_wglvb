% Reading ECCO variables
function var_rsh = ecco_reshape3d(var)

    var_rsh                    = nan(90*3,90*4);
    dimdep                     = size(var,4);
    
    var_rsh(1:90,181:270)      = var(:,:,8);
    var_rsh(91:180,181:270)    = var(:,:,9);
    var_rsh(181:end,181:270)   = var(:,:,10);
    
    var_rsh(1:90,271:end)      = var(:,:,11);
    var_rsh(91:180,271:end)    = var(:,:,12);
    var_rsh(181:end,271:end)   = var(:,:,13);

    for k = 1:dimdep
        var_rsh(1:90,1:90,k)        = flipud(var(:,:,3,k)');
        var_rsh(91:180,1:90,k)      = flipud(var(:,:,2,k)');
        var_rsh(181:end,1:90,k)     = flipud(var(:,:,1,k)');

        var_rsh(1:90,91:180,k)      = flipud(var(:,:,6,k)');
        var_rsh(91:180,91:180,k)    = flipud(var(:,:,5,k)');
        var_rsh(181:end,91:180,k)   = flipud(var(:,:,4,k)');
    end

end