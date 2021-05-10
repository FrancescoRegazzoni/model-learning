function da_test(metamod,obs_err,opt)
    % obs_err [Y*sqrt(T)]
    %% options
    opt.dummy = 0;   
    opt.opt_da.dummy = 0; 
    if ~isfield(opt,'T')
        opt.T = metamod.problem.T;
    end
    if ~isfield(opt,'test')
        opt.test = [];
    end
    if ~isfield(opt.opt_da,'obs_err')
        opt.opt_da.obs_err = obs_err;
    end
        
    %% execution
    obs_err = adapt_dimension(obs_err,metamod.problem.nY);
    if isempty(opt.test)
        opt.test = model_get_random_test(metamod);
    end
    mod = metamodel_particularize(metamod);
    out_HF = model_solve(opt.test,mod);
    out_HF_noised = out_HF;
    out_HF_noised.yy = out_HF.yy + obs_err/sqrt(metamod.dt).*randn(metamod.problem.nY,size(out_HF.yy,2));
    out_da = model_da(metamod,out_HF_noised,opt.opt_da);

    if isfield(mod, 'alpha')
        fprintf('err on alpha:   %1.2e\n',norm(mod.alpha - out_da.xx(metamod.nX+1:end,end)))
    end
    
    figure('units','normalized','outerposition',[0 0 1 1])    
    opt_da_plot.test_clean = out_HF;
    if isfield(mod, 'alpha')
        opt_da_plot.alpha_true = mod.alpha;
    end
    da_plot_output(out_da,out_HF_noised,metamod,opt_da_plot)
%     linewdt = 1.2;
%     figure();
%     nrows = 1;
%     ncols = 1;
%     if metamod.nA > 0, nrows = 2; end
%     subplot(nrows,ncols,1)
%     plot(out_HF.tt,out_HF.yy,'-','linewidth',linewdt); hold on; ax = gca; ax.ColorOrderIndex = 1;
%     plot(out_HF.tt,out_HF_noised.yy,'-','linewidth',1); hold on; ax = gca; ax.ColorOrderIndex = 1;
%     if isequal(metamod.output_type,'insidestate')
%         plot(out_HF.tt,out_da.xx(1:metamod.problem.nY,:),':','linewidth',linewdt); hold on; ax = gca; ax.ColorOrderIndex = 1;
%         plot(out_HF.tt(2:end),out_da.xx_pred(1:metamod.problem.nY,2:end),'--','linewidth',linewdt); hold on; ax = gca; ax.ColorOrderIndex = 1;
%     else
%         plot(out_HF.tt(2:end),out_da.yy_pred(:,2:end),'--','linewidth',linewdt); hold on; ax = gca; ax.ColorOrderIndex = 1;
%     end
%     % axis([0 Tmax 0 3])
%     if metamod.nA > 0
%         subplot(nrows,ncols,2)
% %         plot(out_HF.tt,out_da.xx(metamod.problem.nY+1:end,:),':','linewidth',linewdt); hold on; ax = gca; ax.ColorOrderIndex = 1;
% %         plot(out_HF.tt(2:end),out_da.xx_pred(metamod.problem.nY+1:end,2:end),'--','linewidth',linewdt); hold on; ax = gca; ax.ColorOrderIndex = 1;
%         plot(out_HF.tt,out_da.xx(metamod.nX+1:end,:),':','linewidth',linewdt); hold on; ax = gca; ax.ColorOrderIndex = 1;
%         plot(out_HF.tt(2:end),out_da.xx_pred(metamod.nX+1:end,2:end),'--','linewidth',linewdt); hold on; ax = gca; ax.ColorOrderIndex = 1;
%         plot(out_HF.tt([1 end]),mod.alpha*[1 1])
%     end
%     % axis([0 Tmax a_min a_max])
    
    
end