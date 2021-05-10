function [a0_mean,a0_wide] = metamodel_get_stat_alpha(metamod)

    [a0_mean,a0_wide] = metamodel_get_stat(metamod,'alpha_min','alpha_max','alpha_getrandom',metamod.nA);
    
end