function dataset_learned = model_get_train_fit(ANNmod, opt)

    opt.dummy = 0;
    if ~isfield(opt, 'do_plot')
        opt.do_plot = 1;
    end
    if ~isfield(opt, 'opt_plot')
        opt.opt_plot.dummy = 0;
        opt.opt_plot.size_opt.top_margin = 20;
    end
    
    %%
    d = datasetcouple_get(ANNmod.datasets_def);
    dataset_learned = {};
    for i = 1:length(d.train)
        if ANNmod.problem.samples_variability || ~ANNmod.problem.fixed_x0
            alphas = [];
            IC = [];
            if ANNmod.problem.samples_variability
                alphas = ANNmod.alpha_learned(:,i);
            end
            if ~ANNmod.problem.fixed_x0
                IC = ANNmod.IC_learned(:,i);
            end
            mod_part = metamodel_particularize(ANNmod, IC, alphas);
        else
            mod_part = ANNmod;
        end
        dataset_learned{i} = model_solve(d.train{i}, mod_part, struct('save_x',1));
    end
    
    %%
    if opt.do_plot
        
        if ~isfield(opt.opt_plot, 'plot_x') 
            opt.opt_plot.plot_x = 1;
        end
        if isequal(ANNmod.output_type, 'insidestate')
            if ANNmod.problem.nY == ANNmod.problem.nY
                opt.opt_plot.plot_x = 0;
            end
            opt.opt_plot.plot_x_idxs = ANNmod.problem.nY+1:ANNmod.nX;
        end
        
        dataset_plot(dataset_learned, ANNmod.problem, opt.opt_plot);
    end
    
end