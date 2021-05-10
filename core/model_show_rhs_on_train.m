function model_show_rhs_on_train(mod, opt)

    opt.dummy = 0;
    if ~isfield(opt, 'rhs_limits')
        opt.rhs_limits = [-inf inf];
    end
    if ~isfield(opt, 'size_opt')
        opt.size_opt.dummy = 0;
        opt.size_opt.top_margin = 20;
        opt.size_opt.margin_plot_v = 40;
        opt.size_opt.subplot_width = 300;
        opt.size_opt.subplot_heigth = 200;
    end
    
    %% checks
    if mod.problem.nU ~= 1 || mod.nX ~= 1
        warning('wrong dimensions in model')
        return
    end
    
    %% execution
    nptx = 40;
    nptu = 40;
    u_min = mod.problem.u_min;
    u_max = mod.problem.u_max;
    if isequal(mod.output_type, 'insidestate')
        x_min = mod.problem.y_min;
        x_max = mod.problem.y_max;
    else
        error('not (yet) implemented')
    end
    xx = linspace(x_min,x_max,nptx);
    uu = linspace(u_min,u_max,nptu);
    [UU,XX] = meshgrid(uu,xx);
    FF = zeros(nptu,nptx);  
        
    d = datasetcouple_get(mod.datasets_def);
    nS = length(d.train);
    
    num_rows = ceil(sqrt(nS));
    num_cols = ceil(nS/num_rows);
    win_sizes = plotting_get_axes(num_rows,num_cols, opt.size_opt);
    figure('units','pixel','position',[100 100 win_sizes.width win_sizes.heigth]);
    
    for iS = 1:nS
        if mod.problem.samples_variability || ~problem.fixed_x0
            alphas = [];
            IC = [];
            if mod.problem.samples_variability
                alphas = mod.alpha_learned(:,iS);
            end
            if ~mod.problem.fixed_x0
                IC = mod.IC_learned(:,iS);
            end
            mod_part = metamodel_particularize(mod, IC, alphas);
        else
            mod_part = mod;
        end
                
        for ix = 1:nptx
            for iu = 1:nptu
                FF(ix,iu) = mod_part.f(xx(ix),uu(iu));
            end
        end
        
        out = model_solve(d.train{iS}, mod_part, struct('save_x',1));
    
        subplot('Position',win_sizes.axs{iS}.coord_norm)
        imagesc(uu,xx,FF)
        hold on
        contour(UU,XX,FF,[0 0],'k','linewidth',2)
        caxis(opt.rhs_limits)
        colorbar()
        if isequal(mod.output_type, 'insidestate')
            plot(d.train{iS}.uu,d.train{iS}.yy,'g-','linewidth',1.5)
            plot(d.train{iS}.uu(:,1),d.train{iS}.yy(:,1),'gx','linewidth',1.5)
        end
        plot(out.uu,out.xx,'r-','linewidth',1.5)
        plot(out.uu(:,1),out.xx(:,1),'rx','linewidth',1.5)
        
        if isfield(d.train{iS}, 'label')
            label = d.train{iS}.label;
        else
            label = sprintf('test %d', iS);
        end
        title(label)
        
        set(gca,'YDir','normal')
            
        pause(1e-16)
    end
    
end