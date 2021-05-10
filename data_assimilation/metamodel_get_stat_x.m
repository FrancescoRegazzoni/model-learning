function [x_mean,x_wide] = metamodel_get_stat_x(metamod)

    [x_mean,x_wide] = metamodel_get_stat(metamod,'x_min','x_max','x_getrandom',metamod.nX);
    
end