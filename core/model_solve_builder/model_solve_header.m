function output = model_solve(test,model,opt)
    % MODEL_SOLVE solves a given test with a given model.
    %    test.tt: time (vector, or just extrema)
    %    test.uu: input (vector or handle function)
    %    test.yy: exact output (vector, optional)
    %
    % NB: this function is generated automatically by the script
    % "build_model_solve.m", for efficiency reasons. Do not modify it
    % directly. To modify it:
    %  - enter into folder model_solve_builder;
    %  - modifiy the desired section (model_solve_header, model_solve_footer,
    %    etc.);
    %  - run "build_model_solve.m"
    %  - copy the script "model_solve.m" in this folder into "core"
    
    
    opt.dummy = 0;
    
    if model.problem.metaproblem && ~model.problem.particularized
        error('The associated problem is a non-particularized metaproblem!')
    end
    
    if length(test) > 1
        for iS = 1:length(test)
            output{iS} = model_solve(test{iS},model,opt);
        end
        return
    end
    
    if ~isfield(model,'blackbox')
        model.blackbox = 0;
    end    
    
    %% setting default parameters
    if ~isfield(opt,'save_x')
        opt.save_x = 0;
    end
    if ~isfield(opt,'error_compute')
        opt.error_compute = 0;
    end
    if ~isfield(opt,'save_x_freq')
        opt.save_x_freq = 1;
    end
    if ~isfield(opt,'save_y_ex')
        opt.save_y_ex = 1;
    end
    if ~isfield(opt,'verbose')
        opt.verbose = 0;
    end
    if ~isfield(opt,'do_plot')
        opt.do_plot = 0;
    end
    if ~isfield(opt,'do_plot_x')
        opt.do_plot_x = 0;
    end
    if ~isfield(opt,'interpolation_mode_u')
        opt.interpolation_mode_u = 'pointwise';
%         opt.interpolation_mode_u = 'mean_forward';
    end
    if ~isfield(opt,'interpolation_mode_y')
        opt.interpolation_mode_y = 'pointwise';
    end
    
    
    if model.blackbox
        output = model.forward_function(test,opt);
    else

        %% custom initialization
        switch model.advance_type
            case 'nonlinear_explicit'
                if isfield(model,'u_implicit')
                    if model.u_implicit
                        error('advance_type requires u_implicit = 0')
                    end
                end
                model.u_implicit = 0;
        end
        
        %% building time vector
        tt = test.tt(1):model.dt:test.tt(end);
        dtt = [0 tt(2:end)-tt(1:end-1)];

        if ~isfield(test,'tt_y')
            test.tt_y = test.tt;
        end
        idx_tt_y_first = find(tt >= test.tt_y(1),1,'first');
        idx_tt_y_last  = find(tt <= test.tt_y(end),1,'last');
        tt_y = tt(idx_tt_y_first:idx_tt_y_last);

        %% setting constants
        nT = length(tt);
        nT_y = length(tt_y);
        nX = model.nX;
        nY = model.problem.nY;
        nU = model.problem.nU;

        if strcmp(model.advance_type,'ANN')
            if model.nX == nY
                opt.do_plot_x = 0;
            end
        end

        %% building input/output vectors
        if nU > 0        
            if isa(test.uu,'function_handle')
                ufunc = test.uu;
                uu = ufunc(tt);
            else
                uu = interp_time_series(test.tt,test.uu,tt,struct('mode',opt.interpolation_mode_u));
%                 if isequal(tt,test.tt)
%                     uu = test.uu;
%                 else
%                     for iU = 1:size(test.uu,1)
%                         uu(iU,:) = interp1(test.tt,test.uu(iU,:),tt);
%                     end
%                 end
            end

            if model.u_implicit
                uu_eff = uu;
            else
                uu_eff = [zeros(nU,1) uu(:,1:end-1)];
            end
        else
            % I set it to zero, otherwise I get errors while extracting the current u
            uu_eff = zeros(1,nT);
        end

        if isfield(test,'yy') && (opt.error_compute || opt.do_plot || opt.save_y_ex)
            yy_ex = interp_time_series(test.tt_y,test.yy,tt_y,struct('mode',opt.interpolation_mode_y));
%             if isequal(tt_y,test.tt_y)   
%                 yy_ex = test.yy; 
%             else         
%                 for iY = 1:size(test.yy,1)
%                     yy_ex(iY,:) = interp1(test.tt_y,test.yy(iY,:),tt_y);
%                 end
%             end        
        end

        %% initialiation 
        x = model.x0;

        if opt.save_x || opt.do_plot_x
            xx = x;
        end

        yy = model_get_output(model,x);
        
        %% time loop  
        timeInit = tic();