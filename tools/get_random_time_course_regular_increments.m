function u = get_random_time_course_regular_increments(tt,opt)
    % Returns a randomic time course defined on the time 'tt', with the 
    % regular_increments algorithm.

    if ~isfield(opt,'u0')
        %opt.u0 = 0;
        opt.u0 = opt.umin + (opt.umax-opt.umin)*rand();
    end
    if ~isfield(opt,'regDegree')
        opt.regDegree = 2;
    end
    if ~isfield(opt,'fading_coeff')
        opt.fading_coeff = ones(opt.regDegree-1,1);
    end
    if ~isfield(opt,'movingmeansize')
        opt.movingmeansize = 3;
    end

    nT = length(tt);
    u = zeros(1,nT);

    scales = (opt.umax-opt.umin)/2 ./ opt.time_scale.^(0:opt.regDegree-1)';

    y = scales.*randn(opt.regDegree,1);
    y(1) = opt.u0;
    u(1) = y(1);
    for iT = 2:nT
        if y(1) <= opt.umin
            y(2) = abs(y(2));
            %y(3:end) = 0;
        elseif y(1) >= opt.umax
            y(2) = -abs(y(2));
            %y(3:end) = 0;
        end
        if ~isempty(opt.dudtmin)
            if y(2) <= opt.dudtmin
                y(3) = abs(y(3));
            end
        end
        if ~isempty(opt.dudtmax)
            if y(2) >= opt.dudtmax
                y(3) = -abs(y(3));
            end      
        end

        for iDer = 1:opt.regDegree-1
            y(iDer) = y(iDer)*opt.fading_coeff(iDer) + (tt(iT)-tt(iT-1))*y(iDer+1);
        end
        y(end) = scales(end) * randn();
        u(iT) = y(1); 
    end

    u = min(opt.umax,max(opt.umin,u));
    %figure();plot(tt,u);pause();close()
    if opt.movingmeansize > 0
        u = movmean(u,opt.movingmeansize);
    end
        
end