function st = steady_state_generate(model,u_vals,opt)
    % Computes the steady states associates with the inputs 'u_vals'. 
    % Optionally, it stores the results in the file opt.outFile.

    %% options
    opt.dummy = 0;
    if ~isfield(opt,'do_save')
        opt.do_save = 0;
    end
    if ~isfield(opt,'continue_x')
        opt.continue_x = 0;
    end
    if ~isfield(opt,'verbose')
        opt.verbose = 0;
    end
    
    %% initialization
    if opt.do_save 
        baseopt = get_base_options();
        filename = [baseopt.BaseDir '/' model.problem.dir_data '/' opt.outFile];
    end
    
    %% dataset generation
    nS = size(u_vals,2);
    y_vals = zeros(model.problem.nY,nS);
    for iS = 1:nS
        if opt.verbose
            fprintf('solving test %d of %d...\n',iS,nS)
        end
        opt_solve.dummy = 0;
        if opt.continue_x && iS>1
            opt_solve.x0 = x;
        end
        [y_vals(:,iS), x] = model_get_steady_state(model,u_vals(:,iS),opt_solve);     
    end
    
    %% saving file
    if opt.do_save
        save(filename,'u_vals','y_vals');
    end
    
    st.u_vals = u_vals;
    st.y_vals = y_vals;
end