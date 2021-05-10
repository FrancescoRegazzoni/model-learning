function da_plot_output(out_da,test,metamod,opt)

    opt.dummy = 0;
    if ~isfield(opt,'alpha_true')
        opt.alpha_true = [];
    end
    if ~isfield(opt,'test_clean')
        opt.test_clean = [];
    end
    if ~isfield(opt,'title')
        opt.title = '';
    end

    linewdt = 1.2;
    nrows = 2;
    ncols = 1;
    if metamod.nA > 0
        nrows = nrows+1; 
    end
    subplot(nrows,ncols,1)
    if ~isempty(opt.test_clean)
        plot(opt.test_clean.tt_y,opt.test_clean.yy,'k-','linewidth',2); hold on; ax = gca; ax.ColorOrderIndex = 1;
    end
    if ~isfield(test,'tt_y')
        test.tt_y = test.tt;
    end
    plot(test.tt_y,test.yy,'-','linewidth',linewdt); hold on; ax = gca; ax.ColorOrderIndex = 1;
    if isequal(metamod.output_type,'insidestate')
        plot(out_da.tt,out_da.xx(1:metamod.problem.nY,:),':','linewidth',linewdt); hold on; ax = gca; ax.ColorOrderIndex = 1;
        plot(out_da.tt(2:end),out_da.xx_pred(1:metamod.problem.nY,2:end),'--','linewidth',linewdt); hold on; ax = gca; ax.ColorOrderIndex = 1;
        
        if isfield(out_da,'PP')
            band = 3*sqrt(out_da.PP(1:metamod.problem.nY,:));
            up = out_da.xx(1:metamod.problem.nY,:) + band;
            dw = out_da.xx(1:metamod.problem.nY,:) - band;
            cmap = get(0, 'DefaultAxesColorOrder');
            for iX=1:metamod.problem.nY
                col = cmap(1 + mod(iX-1,size(cmap,1)),:);
                patch([out_da.tt out_da.tt(end:-1:1)],[dw(iX,:) up(iX,end:-1:1)], col, 'EdgeColor','none'); hold on
                alpha(.3) 
            end
        end
        
    else
        plot(out_da.tt(2:end),out_da.yy_pred(:,2:end),'--','linewidth',linewdt); hold on; ax = gca; ax.ColorOrderIndex = 1;
    end
    if ~isempty(opt.title)
        title(opt.title)
    elseif isfield(test,'label')
        title(test.label)
    end
    hold off
    
    subplot(nrows,ncols,2)    
    plot(out_da.tt,out_da.xx(1:metamod.nX,:),':','linewidth',linewdt); hold on; ax = gca; ax.ColorOrderIndex = 1;
    plot(out_da.tt(2:end),out_da.xx_pred(1:metamod.nX,2:end),'--','linewidth',linewdt); hold on; ax = gca; ax.ColorOrderIndex = 1;
    if isfield(out_da,'PP')
        band = 3*sqrt(out_da.PP(1:metamod.nX,:));
        up = out_da.xx(1:metamod.nX,:) + band;
        dw = out_da.xx(1:metamod.nX,:) - band;
        cmap = get(0, 'DefaultAxesColorOrder');
        for iX=1:metamod.nX
            col = cmap(1 + mod(iX-1,size(cmap,1)),:);
            patch([out_da.tt out_da.tt(end:-1:1)],[dw(iX,:) up(iX,end:-1:1)], col, 'EdgeColor','none'); hold on
            alpha(.3) 
        end
    end
    hold off
    
    % axis([0 Tmax 0 3])
    if metamod.nA > 0
        subplot(nrows,ncols,3)
    %         plot(out_HF.tt,out_da.xx(metamod.problem.nY+1:end,:),':','linewidth',linewdt); hold on; ax = gca; ax.ColorOrderIndex = 1;
    %         plot(out_HF.tt(2:end),out_da.xx_pred(metamod.problem.nY+1:end,2:end),'--','linewidth',linewdt); hold on; ax = gca; ax.ColorOrderIndex = 1;
        leg_lines = plot(out_da.tt,out_da.xx(metamod.nX+1:end,:),':','linewidth',linewdt); hold on; ax = gca; ax.ColorOrderIndex = 1;
        plot(out_da.tt(2:end),out_da.xx_pred(metamod.nX+1:end,2:end),'--','linewidth',linewdt); hold on; ax = gca; ax.ColorOrderIndex = 1;
        if ~isempty(opt.alpha_true)
             leg_lines = plot(out_da.tt([1 end]),(opt.alpha_true*[1 1])');
        end
        
        if isfield(out_da,'PP')
            band = 3*sqrt(out_da.PP(metamod.nX+1:end,:));
            up = out_da.xx(metamod.nX+1:end,:) + band;
            dw = out_da.xx(metamod.nX+1:end,:) - band;
            cmap = get(0, 'DefaultAxesColorOrder');
            for iX=1:metamod.nA
                col = cmap(1 + mod(iX-1,size(cmap,1)),:);
                patch([out_da.tt out_da.tt(end:-1:1)],[dw(iX,:) up(iX,end:-1:1)], col, 'EdgeColor','none'); hold on
                alpha(.3) 
            end
        end
    
        if isfield(metamod,'params_names')
            if ~isempty(metamod.params_names)
                legend(leg_lines, metamod.params_names,'Location','northwest','Orientation','horizontal')
    %             legend boxoff
            end
        end
        hold off
    end
    % axis([0 Tmax a_min a_max])
%     pause(1e-16)