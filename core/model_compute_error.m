function out_err = model_compute_error(model,tests,opt)
    % Computes the error of 'model' on the dataset 'tests'.
    
    opt.dummy = 0;
    if ~isfield(opt,'verbose')
        opt.verbose = 1;
    end
    if ~isfield(opt,'verbose_eachtest')
        opt.verbose_eachtest = 0;
    end
    if ~isfield(opt,'do_plot')
        opt.do_plot = 0;
    end
    if ~isfield(opt,'do_plot_x')
        opt.do_plot_x = 0;
    end
    if ~isfield(opt,'pause_eachtest')
        opt.pause_eachtest = opt.do_plot;
    end
    
    opt_solve.error_compute = 1;
    opt_solve.verbose = opt.verbose_eachtest;
    opt_solve.do_plot_x = opt.do_plot_x;
    
    if opt.do_plot
        opt_solve.do_plot = 1;
        figure()
    end
    
    err_L2_vec = [];
    norm_L2_vec = [];
    times = [];
    
    if ischar(tests)       
        dataset_def.problem = model.problem;
        dataset_def.type = 'file';
        dataset_def.source = tests;
        tests = dataset_get(dataset_def); 
    end
    
    for iTest=1:length(tests)
        output = model_solve(tests{iTest},model,opt_solve);
        err_L2_vec = [err_L2_vec output.err_L2];
%         norm_L2_vec = [norm_L2_vec sqrt(mean(output.yy(:).^2))];
        norm_L2_vec = [norm_L2_vec get_norm_L2_time(output.tt,output.yy_ex - model.problem.y0_norm_computation)];
        times = [times output.time_norm];
        if opt.pause_eachtest && iTest<length(tests)
            pause()
        end
    end
    
    out_err.err_dataset_L2 = sqrt(mean(err_L2_vec.^2));
    out_err.norm_dataset_L2 = sqrt(mean(norm_L2_vec.^2));
    
    out_err.err_dataset_L2_norm = out_err.err_dataset_L2/out_err.norm_dataset_L2;
    out_err.time_mean = mean(times);
        
    if opt.verbose
        fprintf('overall err L2: %1.2e\n',out_err.err_dataset_L2_norm)
        fprintf('mean time: %1.2e s\n',out_err.time_mean)
    end
    
end