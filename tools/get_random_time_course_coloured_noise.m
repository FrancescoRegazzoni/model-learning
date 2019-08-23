function xx = get_random_time_course_coloured_noise(tt,opt)
    % Returns a randomic time course defined on the time 'tt', with the 
    % coloured_noise algorithm.

    if ~isfield(opt,'exponent')
        opt.exponent = 2;
    end
    if ~isfield(opt,'reduction_range')
        opt.reduction_range = 10;
    end
    
    nT = length(tt);

    corrfunc = @(t) (opt.umax-opt.umin)/opt.reduction_range*exp(-(t/opt.time_scale).^opt.exponent);

    A = eye(nT)*corrfunc(0)^2;
    for i=1:nT-1
        dd = corrfunc(tt(1+i:end)-tt(1:end-i)).^2;
        A = A + diag(dd,i) + diag(dd,-i);
    end
    S = sqrt(A);

    found = 0;
    while ~found
        zz = randn(1,nT);
        xx = .5*(opt.umin+opt.umax) + (S*zz')';
        vv = (xx(2:end)-xx(1:end-1))./(tt(2:end)-tt(1:end-1));
        if max(xx) <= opt.umax && min(xx) >= opt.umin
            found = 1;
        end
    end

end