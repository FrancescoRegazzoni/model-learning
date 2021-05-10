function ret = da_estimate_predict(problem,model_learned,opt)
    %% Parameters
    opt.dummy = 0;
    if ~isfield(opt,'dataset')
        opt.dataset = '';
    end
    if ~isfield(opt,'dataset_add_noise')
        opt.dataset_add_noise = 0;
    end
    if ~isfield(opt,'mod_HF')
        opt.mod_HF = [];
    end
    if ~isfield(opt,'use_training_dataset')
        opt.use_training_dataset = 0;
    end
    if ~isfield(opt,'pause_each_test')
        opt.pause_each_test = 0;
    end
    if ~isfield(opt,'do_plot')
        opt.do_plot = 0;
    end
    if ~isfield(opt,'do_plot_error')
        opt.do_plot_error = 1;
    end
    if ~isfield(opt,'do_plot_da')
        opt.do_plot_da = 0;
    end
    if ~isfield(opt,'do_models_comparison')
        opt.do_models_comparison = 0;
    end
    if ~isfield(opt,'T_tot')
        opt.T_tot = problem.T;
    end
    if ~isfield(opt,'fraction_obs')
        opt.fraction_obs = .5;
    end
    if ~isfield(opt,'use_alpha_to_alpha')
        opt.use_alpha_to_alpha = 1;
    end
    if ~isfield(opt,'obs_err') % [Y*sqrt(T)]
        % nY times num_cases (num_cases is the number of error cases 
        % to test. if greater than one, a comparison is performed).
        opt.obs_err = 1e-1;
    end
    if ~isfield(opt,'da_obs_err') % [Y*sqrt(T)]
        opt.da_obs_err = opt.obs_err;
    end
    if ~isfield(opt,'da_mod_err') % [Y/sqrt(T)]
        opt.da_mod_err = 0;
    end
    if ~isfield(opt,'reference_qty') %obs_err, da_obs_err, da_mod_err, noise_train
        opt.reference_qty = 'obs_err';
    end
    if ~isfield(opt,'n_tests')
        opt.n_tests = 1e3;
    end
    if ~isfield(opt,'dataset_alpha')
        opt.dataset_alpha = [];
    end
    if ~isfield(opt,'out_alpha')
        opt.out_alpha = 0;
    end
    if ~isfield(opt,'out_last_data_assimilation_prediction')
        opt.out_last_data_assimilation_prediction = 0;
    end
    
    %% Check
    if isempty(opt.mod_HF) && isempty(opt.dataset) && ~opt.use_training_dataset
        error('At least one among HF_mod, dataset and use_training_dataset bust be passed')
    end
    if ~isempty(opt.dataset_alpha) && isempty(opt.dataset)
        error('dataset_alpha can be used only in combination with a dataset')
    end
    
    if opt.use_training_dataset
        opt.dataset = model_learned.datasets_def.source_train;
        opt.mod_HF = model_learned;
        if problem.samples_variability
            opt.dataset_alpha = model_learned.alpha_learned;
            model_learned.alpha_original = model_learned.alpha_learned;
            model_learned.alpha_to_alpha = @(x) x;
        end
    end
    
    %% Initialization    
    if istable(model_learned)
        num_models = size(model_learned,1);
        fprintf('reading models...\n')
        for i_mod = 1:num_models
            model_list{i_mod} = read_model_fromfile(problem,model_learned.dir{i_mod});
            noise_list(1,i_mod) = model_learned.noise(i_mod);
        end
        % model used for inizialization (assuming that all the models share the same properties)
        model_learned = model_list{1}; 
    else
        num_models = 1;
    end
    
    if problem.samples_variability
        nA = model_learned.nA;
    else
        nA = 0;
    end
    nX = model_learned.nX;
    nY = model_learned.problem.nY;
    
    if ~isempty(opt.dataset)
        if ~ opt.dataset_add_noise
            opt.obs_err = 0;
        end
        if ischar(opt.dataset)
            samples = dataset_get(struct('problem',problem,'type','file','source',opt.dataset));
        else
            samples = opt.dataset;
        end
        n_tests = length(samples);
        if ~isempty(opt.dataset_alpha) 
            if size(opt.dataset_alpha,2) ~= n_tests
                error('the size of dataset_alpha must be equal to the dataset size')
            end
        end
    else
        n_tests = opt.n_tests;
    end
    
    n_cases = max([size(opt.obs_err,2),size(opt.da_obs_err,2),size(opt.da_mod_err,2),num_models]);
    opt.obs_err    = opt.obs_err   .*ones(1,n_cases);
    opt.da_obs_err = opt.da_obs_err.*ones(1,n_cases);
    opt.da_mod_err = opt.da_mod_err.*ones(1,n_cases);
    
    err_x_rel_cum = zeros(1,n_cases);
    err_alpha_rel_cum = zeros(1,n_cases);
    err_pred_rel_cum = zeros(1,n_cases);
    
    T_estimate = opt.fraction_obs*opt.T_tot;

    alpha_to_alpha_available = 1;
    if model_learned.nA > 0
        if isempty(opt.mod_HF)
            auto_prediction = 0;
        else
            auto_prediction = isequal(opt.mod_HF,model_learned);
        end
        if auto_prediction
            alpha_to_alpha = @(a) a;
        elseif opt.use_alpha_to_alpha && isfield(model_learned,'alpha_to_alpha')
            alpha_to_alpha = model_learned.alpha_to_alpha;   
        elseif opt.use_alpha_to_alpha && isfield(model_learned,'get_alpha_to_alpha_handler')
            alpha_to_alpha = model_learned.get_alpha_to_alpha_handler(opt.mod_HF); 
        else
            alpha_to_alpha_available = 0;
        end
    else
        alpha_to_alpha_available = 0;
    end
        

    if opt.do_plot
        fig_da = figure();
        if problem.samples_variability
            if isfield(model_learned,'alpha_min') && isfield(model_learned,'alpha_max')
                alpha_min_fig = min(model_learned.alpha_min);
                alpha_max_fig = max(model_learned.alpha_max);
            else
                [a0_mean,a0_wide] = metamodel_get_stat_alpha(model_learned);
                alpha_min_fig = min(a0_mean - 3*a0_wide);
                alpha_max_fig = max(a0_mean + 3*a0_wide);
            end
            alpha_min_fig = alpha_min_fig - .1*(alpha_max_fig-alpha_min_fig);
            alpha_max_fig = alpha_max_fig + .1*(alpha_max_fig-alpha_min_fig);
            if nA == 1
                if auto_prediction
                    fig_alpha = figure();
                    [a0_mean,a0_wide] = metamodel_get_stat_alpha(opt.mod_HF);
                    aa_plot = linspace(a0_mean-a0_wide*sqrt(12)/2,a0_mean+a0_wide*sqrt(12)/2,1e2);
                    plot(aa_plot,aa_plot,'k-')
                    hold on
                elseif isfield(model_learned, 'alpha_original')
                    fig_alpha = model_alpha_plot(model_learned);
                end
            elseif nA == 2
                if isfield(model_learned, 'alpha_original')
                    fig_alpha = model_alpha_plot(model_learned);
                else
                    fig_alpha = figure();
                end
            end
        end
% %         x0_norm = ones(model_learned.nX,1); %TODO
%         [x0_norm, ~] = metamodel_get_stat_x0(model_learned);
% %         x_min_fig = max(mod_HF.problem.y0_min ./ x0_norm);
% %         x_max_fig = min(mod_HF.problem.y0_max ./ x0_norm);
%         x_min_fig = 0; % TODO
%         x_max_fig = 3; % TODO
        y_norm = .5*(model_learned.problem.y_max-model_learned.problem.y_min);
        y_min_fig = min(model_learned.problem.y_min ./ y_norm);
        y_max_fig = max(model_learned.problem.y_max ./ y_norm);
        
        if opt.do_models_comparison
            fig_mode_comparison = figure('units','normalized','outerposition',[0 0 1 1]);
        end
    end
        
    if n_cases > 1 && opt.do_plot_error
        if alpha_to_alpha_available
            fig_comparison = figure('units','pixel','outerposition',[200 200 1200 400]);
            nr_cols = 3;
        else
            fig_comparison = figure('units','pixel','outerposition',[200 200 800 400]);
            nr_cols = 2;
        end        
        switch opt.reference_qty
            case 'obs_err' 
                x_values = mean(opt.obs_err,1);
                x_label = 'applied observation error';
            case 'da_obs_err' 
                x_values = mean(opt.da_obs_err,1);
                x_label = 'observation error (in da)';
            case 'da_mod_err' 
                x_values = mean(opt.da_mod_err,1);
                x_label = 'model error (in da)';
            case 'noise_train' 
                x_values = noise_list;
                x_label = 'train noise';
        end
    end
    
    num_char = '10';
%     da_estimate_predict_printheader();

    if opt.out_alpha
        alpha_estimated = zeros(model_learned.nA,n_tests,n_cases);
        alpha_original = zeros(opt.mod_HF.nA,n_tests,n_cases);
    end
    
    for i_test = 1:n_tests
        for i_case = 1:n_cases
            %% Test generation
            if isempty(opt.dataset)
                test_solve.tt = [0 opt.T_tot];
                mod_HF_part = metamodel_particularize(opt.mod_HF);
                alpha_exact = mod_HF_part.alpha;
                alpha_exact_available = 1;
                out_HF = model_solve(test_solve,mod_HF_part);
            else
                if ~isempty(opt.dataset_alpha) 
                    alpha_exact = opt.dataset_alpha(:,i_test);
                    alpha_exact_available = 1;
                elseif isfield(samples{i_test}, 'alpha')
                    alpha_exact = samples{i_test}.alpha;
                    alpha_exact_available = 1;
                else
                    alpha_exact_available = 0;
                end
                out_HF = samples{i_test};
                if isfield(samples{i_test}, 'yy_ex')
                    out_HF.yy = samples{i_test}.yy_ex;
                end
            end
            
            %% Adding noise
            iTda = max(find(out_HF.tt < T_estimate));
            out_HF_noised.tt = out_HF.tt(:,1:iTda);
            if isfield(out_HF, 'uu')
                out_HF_noised.uu = out_HF.uu(:,1:iTda);
            end
            if isempty(opt.dataset) || opt.dataset_add_noise
        %         out_HF_noised = out_HF;
                out_HF_noised.yy = out_HF.yy(:,1:iTda) + opt.obs_err(:,i_case)/sqrt(opt.mod_HF.dt).*randn(nY,iTda);
            else
                out_HF_noised.yy = samples{i_test}.yy(:,1:iTda);
            end

            %% Data assimilation
%             out_da = model_da(model_learned,out_HF_noised,struct('obs_err',obs_err(i_case)));
            if num_models > 1
                model_learned = model_list{i_case};
            end
            out_da = model_da(model_learned,out_HF_noised,struct('obs_err',opt.da_obs_err(:,i_case),'mod_err',opt.da_mod_err(:,i_case)));
            x_star = out_da.xx(1:nX,end);
            if nA > 0
                alpha_star = out_da.xx(nX+1:end,end);
                if opt.out_alpha
                    alpha_estimated(:,i_test,i_case) = alpha_star;
                    if alpha_exact_available
                        alpha_original(:,i_test,i_case) = alpha_exact;
                    end
                end
            else
                alpha_star = [];
            end
            
            %% Plot (estimation)
            if opt.do_plot
                linewdt = 1.2;
                
                n_rows = 1;
                if nX > nY || ~isequal(model_learned.output_type, 'insidestate')
                    n_rows = n_rows + 1;
                end
                if nA > 0
                    n_rows = n_rows + 1;
                end
                
                band = 3*sqrt(out_da.PP);
                up = out_da.xx + band;
                dw = out_da.xx - band;
            
                figure(fig_da);
                
                idx_row = 1;
                subplot(n_rows,1,idx_row)  
                hold off;
                
                plot(out_HF.tt,out_HF.yy./y_norm,'-','linewidth',linewdt); hold on;
                if max(opt.obs_err(:,i_case)) > 0        
                    ax = gca;
                    ax.ColorOrderIndex = 1;
                    plot(out_HF.tt(1:iTda),out_HF_noised.yy./y_norm,'-','linewidth',linewdt); hold on;
                end
                plot(out_da.tt,out_da.yy./y_norm,'k--','linewidth',linewdt); hold on;
                
%                 if isequal(model_learned.output_type, 'insidestate')
%                     plot(out_da.tt,out_da.xx(1:nY,:)./y_norm,'k--','linewidth',linewdt); hold on;
%                     plot(out_da.tt,out_da.xx(nY+1:nX,:), '--','linewidth',linewdt); hold on;
%                 else
%                     plot(out_da.tt,out_da.xx(1:nX,:),'--','linewidth',linewdt); hold on;
%                 end
                if isequal(model_learned.output_type, 'insidestate')
                    cmap = get(0, 'DefaultAxesColorOrder');
                    for iX=1:nY
                        col = cmap(1 + mod(iX-1,size(cmap,1)),:);
                        patch([out_da.tt out_da.tt(end:-1:1)],[dw(iX,:) up(iX,end:-1:1)]/y_norm(iX), col, 'EdgeColor','none'); hold on
                        alpha(.3) 
                    end
                end                
                axis([0 max(out_HF.tt) y_min_fig y_max_fig])
                ylabel('y')
                xline(T_estimate);
                
                if nX > nY || ~isequal(model_learned.output_type, 'insidestate')                    
                    idx_row = idx_row + 1;
                    subplot(n_rows,1,idx_row)         
                    hold off;
                    if isequal(model_learned.output_type, 'insidestate')
                        offset_idx_x = nY;
                    else
                        offset_idx_x = 0;
                    end
                    
                    plot(out_da.tt,out_da.xx(offset_idx_x+1:nX,:),'k--','linewidth',linewdt); hold on;
                    
                    cmap = get(0, 'DefaultAxesColorOrder');
                    for iX=1:(nX - offset_idx_x)
                        col = cmap(1 + mod(iX-1,size(cmap,1)),:);
                        patch([out_da.tt out_da.tt(end:-1:1)],[dw(offset_idx_x + iX,:) up(offset_idx_x + iX,end:-1:1)], col, 'EdgeColor','none'); hold on
                        alpha(.3) 
                    end
                    
                    axis([0 max(out_HF.tt) -inf inf])
                    ylabel('x')
                    xline(T_estimate);
                end
                                
                if nA > 0                    
                    idx_row = idx_row + 1;
                    subplot(n_rows,1,idx_row)    
                    hold off;
                    if alpha_to_alpha_available && alpha_exact_available
                        plot(out_HF.tt(),alpha_to_alpha(alpha_exact)*ones(1,length(out_HF.tt)),'--','linewidth',.8); hold on; 
                        ax = gca; ax.ColorOrderIndex = 1;
                    end
                    plot(out_da.tt,out_da.xx(nX+1:end,:),'-','linewidth',linewdt); hold on;
                    cmap = get(0, 'DefaultAxesColorOrder');
                    for iA=1:model_learned.nA
                        col = cmap(1 + mod(iA-1,size(cmap,1)),:);
                        patch([out_da.tt out_da.tt(end:-1:1)],[dw(nX+iA,:) up(nX+iA,end:-1:1)], col, 'EdgeColor','none'); hold on
                        alpha(.3) 
                    end
                
                    axis([0 max(out_HF.tt) alpha_min_fig alpha_max_fig])
                    ylabel('\alpha')
                    xline(T_estimate);

                    if alpha_exact_available
                        if nA == 1
                            figure(fig_alpha)
                            plot(alpha_star,alpha_exact,'k+','MarkerSize',20);
                        elseif nA == 2
                            figure(fig_alpha)
                            subplot(1,2,1)
                            plot(alpha_star(1),alpha_star(2),'k+','MarkerSize',20);
                            text(alpha_star(1),alpha_star(2),sprintf('#%d', i_test))
                            subplot(1,2,2)
                            plot(alpha_exact(1),alpha_exact(2),'k+','MarkerSize',20);
                            text(alpha_exact(1),alpha_exact(2),sprintf('#%d', i_test))
                        end   
                        if opt.do_models_comparison
                            figure(fig_mode_comparison);
                            models_comparison(metamodel_particularize(model_learned,[],alpha_star ), ...
                                              metamodel_particularize(opt.mod_HF,   [],alpha_exact));
                        end 
                    end
                end
            end

            %% prediction
            model_learned_part = metamodel_particularize(model_learned,x_star,alpha_star);
            test_predict.tt = out_HF.tt(:,iTda:end);
            test_predict.tt_y = test_predict.tt;
            test_predict.yy = out_HF.yy(:,iTda:end);
            if isfield(out_HF, 'uu')
                test_predict.uu = out_HF.uu(:,iTda:end);
            end
            out_predict = model_solve(test_predict,model_learned_part, struct('save_x',1));

            %% Errors computation
            % Error on state
            err_x = out_HF.yy(:,iTda) - x_star;
            err_x_rel = norm(err_x)/norm(x_star);
            err_x_rel_cum(i_case) = err_x_rel_cum(i_case) + err_x_rel;
            err_x_rel_list(i_test,i_case) = err_x_rel;
            % Error on alpha
            if alpha_to_alpha_available && alpha_exact_available
                err_alpha = alpha_to_alpha(alpha_exact) - alpha_star;
                err_alpha_rel = norm(err_alpha)/norm(alpha_star);
                err_alpha_rel_cum(i_case) = err_alpha_rel_cum(i_case) + err_alpha_rel;
                err_alpha_rel_list(i_test,i_case) = err_alpha_rel;
            end
            % Error on prediction
            yy_ex = out_HF.yy(:,iTda:end);
            yy_es = out_predict.yy;
            if size(yy_ex,2) < size(yy_es,2)
                yy_es = interp_time_series(out_predict.tt,yy_es,out_HF.tt(:,iTda:end));
            elseif size(yy_ex,2) > size(yy_es,2)
                yy_ex = interp_time_series(out_HF.tt(:,iTda:end),yy_ex,out_predict.tt);
            end
            err_pred = yy_es - yy_ex;
            err_pred_rel = norm(err_pred)/norm(yy_ex);            
            err_pred_rel_cum(i_case) = err_pred_rel_cum(i_case) + err_pred_rel;
            err_pred_rel_list(i_test,i_case) = err_pred_rel;

            %% Plot (prediction)
            if opt.do_plot
                figure(fig_da)
                
                idx_row = 1;
                subplot(n_rows,1,idx_row)    
                hold on
                plot(out_predict.tt,out_predict.yy./y_norm,'k--','linewidth',2)
                                
                if nX > nY || ~isequal(model_learned.output_type, 'insidestate')                    
                    idx_row = idx_row + 1;
                    subplot(n_rows,1,idx_row)            
                    hold on
                    plot(out_predict.tt, out_predict.xx(offset_idx_x+1:end,:),'k--','linewidth',2)
                end
            end
            
%             fprintf('num tests: %d  obs_err = %1.2e   err = %1.5f    mean err = %1.2e\n',i_test,obs_err(i_case),err_pred_rel,errTot(i_case)/i_test)
%             fprintf('num tests: %d  obs_err = %1.2e   err = %1.5f    mean err = %1.2e\n',i_test,obs_err(i_case),err_pred_rel,errTot(i_case)/i_test)
            
            da_estimate_predict_printheader();
            fprintf(['%' num_char 'd   %' num_char 'g   %' num_char 'f   %' num_char 'f   %' num_char 'f   %' num_char 'f'], i_test,mean(opt.obs_err(i_case)),err_pred_rel,err_pred_rel_cum(i_case)/i_test,err_x_rel,err_x_rel_cum(i_case)/i_test);
            if alpha_to_alpha_available
                fprintf(['  %' num_char 'f   %' num_char 'f'], err_alpha_rel,err_alpha_rel_cum(i_case)/i_test);
            end
            fprintf('\n');
            
            if opt.do_plot
                if opt.pause_each_test
                    pause(); 
                else
                    pause(1e-16);
                end
            end
        end
%         if n_cases > 1 && mod(i_test,10) == 1
        if n_cases > 1 && opt.do_plot_error
            figure(fig_comparison)
%             loglog(obs_err,obs_err,'k--')
%             hold on
%             loglog(obs_err,err_pred_rel_list,'k+')
%             p_max = loglog(obs_err,max(err_pred_rel_list,[],1),'bo-','linewidth',1.2);
%             p_min = loglog(obs_err,min(err_pred_rel_list,[],1),'bo-','linewidth',1.2);
% %             loglog(obs_err,errTot/i_test,'ro-','linewidth',1.2)
%             p_arit_mean = loglog(obs_err,mean(err_pred_rel_list,1),'ro-','linewidth',1.2);
%             p_geom_mean = loglog(obs_err,exp(mean(log(err_pred_rel_list),1)),'go-','linewidth',1.2);
%             grid on
%             xlabel('observation error')
%             ylabel('prediction error')
%             legend([p_arit_mean p_geom_mean p_max],{'arithmetic mean', 'geometric mean', 'min/max'},'Location','southoutside','Orientation','horizontal')
%             axis equal
%             hold off
%             pause(1e-16);
            subplot(1,nr_cols,1)
            da_estimate_predict_makeplot(err_pred_rel_list,'prediction error')
            subplot(1,nr_cols,2)
            da_estimate_predict_makeplot(err_x_rel_list,'state error')
            if alpha_to_alpha_available
                subplot(1,nr_cols,3)
                da_estimate_predict_makeplot(err_alpha_rel_list,'alpha error')
            end
            pause(1e-16);
        end
    end
    
    ret.dummy = 0;
    if n_cases > 1        
        ret.obs_err = opt.obs_err;
        ret.da_obs_err = opt.da_obs_err;
        ret.da_mod_err = opt.da_mod_err;
        if istable(model_learned)
            ret.noise_list = noise_list;
        end
        ret.err_pred_rel_list = err_pred_rel_list;
        ret.err_x_rel_list = err_x_rel_list;
        if alpha_to_alpha_available
            ret.err_alpha_rel_list = err_alpha_rel_list;
        end
    end
    if opt.out_alpha
        ret.alpha_estimated = alpha_estimated;
        ret.alpha_original = alpha_original;
    end
    if opt.out_last_data_assimilation_prediction
        ret.out_da = out_da;
        ret.out_predict = out_predict;
    end

    function da_estimate_predict_printheader()
        fprintf(['%' num_char 's   %' num_char 's   %' num_char 's   %' num_char 's   %' num_char 's   %' num_char 's'], 'tests num.', 'obs. err.', 'pred. err', '(mean)', 'state err.', '(mean)');
        if alpha_to_alpha_available
            fprintf(['  %' num_char 's   %' num_char 's'], 'alpha err.', '(mean)');
        end
        fprintf('\n');
%         sep_str = '----------';
%         fprintf(['%' num_char 's   %' num_char 's   %' num_char 's   %' num_char 's  %' num_char 's  %' num_char 's'], sep_str, sep_str, sep_str, sep_str, sep_str, sep_str);
%         if alpha_to_alpha_available
%             fprintf(['  %' num_char 's  %' num_char 's'], sep_str, sep_str);
%         end
%         fprintf('\n');
    end

    function da_estimate_predict_makeplot(err_list,tit)
%         loglog(x_values',mean(obs_err,1)','k--'); hold on        
        loglog(x_values',err_list','k+'); hold on
        p_max = loglog(x_values',max(err_list,[],1)','bo-','linewidth',1.2);
        p_min = loglog(x_values',min(err_list,[],1)','bo-','linewidth',1.2);
        p_arit_mean = loglog(x_values',mean(err_list,1)','ro-','linewidth',1.2);
        p_geom_mean = loglog(x_values',exp(mean(log(err_list),1))','go-','linewidth',1.2);
        grid on
        xlabel(x_label)
        title(tit)
% %         legend([p_arit_mean p_geom_mean p_max],{'a. mean', 'g. mean', 'min/max'},'Location','southoutside','Orientation','horizontal')
%         legend([p_arit_mean p_geom_mean p_max],{'a.m.', 'g.m.', 'min/max'},'Location','south','Orientation','horizontal')
%         legend boxoff
% %         axis equal
        hold off
    end
end
