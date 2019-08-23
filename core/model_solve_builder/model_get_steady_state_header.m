function [y,x] = model_get_steady_state(model,u_steady,opt)
    
    % NB: this function is generated automatically by the script
    % "build_model_get_steady_state.m", for efficiency reasons. Do not 
    % modify it directly. To modify it:
    %  - enter into folder model_solve_builder;
    %  - modifiy the desired section (model_get_steady_state_header, 
    %    model_get_steady_state_loop_footer,etc.);
    %  - run "build_model_get_steady_state.m"
    %  - copy the script "model_get_steady_state.m" in this folder into "core"

    %% options
    opt.dummy = 0;
    if ~isfield(opt,'eps_x')
        opt.eps_x = 1e-4;
    end
    if ~isfield(opt,'eps_y')
        opt.eps_y = 1e-4;
    end
    if ~isfield(opt,'x0')
        opt.x0 = [];
    end
    
    %% initialization
    if isempty(opt.x0)
        x = model.x0;
    else
        x = opt.x0;
    end
    
    dt = model.dt;
    
    switch model.output_type
        case 'linear' 
            g = @(x) model.G*x+model.g0;
        case 'nonlinear'
            g = model.g;
        case 'insidestate'
            g = @(x) x(1:model.problem.nY);
        otherwise
            error('unknown output type %s',model.output_type)
    end
    
    y = g(x);
    
    y_norm = model.problem.y_max-model.problem.y_min;
    
    %% steady state computation
    steady_state = 0;