function DF = numerical_diff(F,x0,opt)

    if ~isfield(opt,'ep_rel')
        ep_rel = 1e-6;
    end
    if ~isfield(opt,'ep_min')
        ep_min = 1e-20;
    end
    if isfield(opt,'F0')
        F0 = opt.F0;
    else
        F0 = F(x0);
    end
    
    num_x = length(x0);
    num_F = length(F0);
    DF = nan(num_F,num_x);
    for i_x = 1:num_x
        x_new = x0;        
        ep = max(ep_min,abs(x0(i_x))*ep_rel);
        x_new(i_x) = x_new(i_x) + ep;
        DF(:,i_x) = (F(x_new)-F0)/ep;
    end
end