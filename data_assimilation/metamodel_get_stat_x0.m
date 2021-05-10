function [x0_mean,x0_wide] = metamodel_get_stat_x0(metamod)

    if metamod.problem.fixed_x0
        x0_mean = metamod.x0;
        x0_wide = zeros(metamod.nX,1);
    else
        [x0_mean,x0_wide] = metamodel_get_stat(metamod,'x0_min','x0_max','x0_getrandom',metamod.nX);
    end
    
end