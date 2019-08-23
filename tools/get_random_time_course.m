function u = get_random_time_course(tt,opt)
    % Returns a randomic time course defined on the time 'tt'.

    if nargin == 1
        opt.dummy = 0;
    end

    if ~isfield(opt,'dim')
        opt.dim = 1;
    end

    if opt.dim > 1
        u = [];
        optScalar = opt;
        optScalar.dim = 1;
        optScalar.do_plot = 0;
        for id = 1:opt.dim
            if isfield(opt,'umin')
                if length(opt.umin) > 1
                    optScalar.umin = opt.umin(id);
                end
            end
            if isfield(opt,'umax')
                if length(opt.umax) > 1
                    optScalar.umax = opt.umax(id);
                end
            end
            if isfield(opt,'dudtumin')
                if length(opt.dudtumin) > 1
                    optScalar.dudtumin = opt.dudtumin(id);
                end
            end
            if isfield(opt,'dudtmax')
                if length(opt.dudtmax) > 1
                    optScalar.dudtmax = opt.dudtmax(id);
                end
            end
            if isfield(opt,'time_scale')
                if length(opt.time_scale) > 1
                    optScalar.time_scale = opt.time_scale(id);
                end
            end
            if isfield(opt,'exponent')
                if length(opt.exponent) > 1
                    optScalar.exponent = opt.exponent(id);
                end
            end
            if isfield(opt,'reduction_range')
                if length(opt.reduction_range) > 1
                    optScalar.reduction_range = opt.reduction_range(id);
                end
            end
            u = [u; get_random_time_course(tt,optScalar)];
        end
    else    
        
        if ~isfield(opt,'method') % regular_increments; coloured_noise; 
            opt.method = 'regular_increments';
        end
        if ~isfield(opt,'umin')
            opt.umin = -1;
        end
        if ~isfield(opt,'umax')
            opt.umax = 1;
        end
        if ~isfield(opt,'dudtmin')
            opt.dudtmin = [];
        end
        if ~isfield(opt,'dudtmax')
            opt.dudtmax = [];
        end
        if ~isfield(opt,'time_scale')
            opt.time_scale = tt(end)/10;
        end
        
        switch opt.method
            case 'regular_increments'
                u = get_random_time_course_regular_increments(tt,opt);
            case 'coloured_noise'
                u = get_random_time_course_coloured_noise(tt,opt);
            otherwise
                error('method %s not exists.', opt.method)
        end
    end


    if isfield(opt,'do_plot')
        if opt.do_plot
            figure();
            subplot(2,1,1)
            plot(tt,u);
            subplot(2,1,2)
            plot(tt(1:end-1),(u(:,2:end)-u(:,1:end-1))./(tt(2:end)-tt(1:end-1)));
            pause(1e-16);
        end
    end

end