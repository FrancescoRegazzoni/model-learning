function dataset_plot(dataset,problem,opt)
    % Plots a dataset.

    fprintf('plotting dataset...')
    %% setting default values
    opt.dummy = 0;
    opt.size_opt.dummy = 0;
    if ~isfield(opt,'figure_size') 
        opt.figure_size = [];
    end
    if ~isfield(opt,'one_sample_per_time') 
        opt.one_sample_per_time = 0;
    end
    if ~isfield(opt,'idx_samples') 
        opt.idx_samples = [];
    end
    if ~isfield(opt,'plot_x') 
        opt.plot_x = 0;
    end
    if ~isfield(opt,'plot_x_idxs') 
        opt.plot_x_idxs = -1;
    end
    if ~isfield(opt,'do_save') 
        opt.do_save = 0;
    end
    if ~isfield(opt,'do_save_external_legend') 
        opt.do_save_external_legend = 0;
    end
    if ~isfield(opt,'filename') 
        opt.filename = 'samples';
    end
    %%%%%%%% num rows/cols
    if ~isfield(opt,'num_cols') 
        opt.num_cols = 0;
    end
    if ~isfield(opt,'num_rows') 
        opt.num_rows = 0;
    end
    %%%%%%%% aspect
    if ~isfield(opt,'automargin')
        opt.automargin = 0;
    end
    if ~isfield(opt,'uy_same_axis')
        opt.uy_same_axis = 0;
    end
    if ~isfield(opt,'plot_legend')
        opt.plot_legend = 0;
    end
    if ~isfield(opt,'show_x_ticks') % 0 = never; 1 = always; -1 = on borders
        opt.show_x_ticks = 1;
    end
    if ~isfield(opt,'show_y_ticks') % 0 = never; 1 = always; -1 = on borders
        opt.show_y_ticks = -1;
    end
    if ~isfield(opt,'y_label') 
        opt.y_label = 'y';
    end
    if ~isfield(opt,'u_lim') 
        if problem.nU > 0
            opt.u_lim = [min(problem.u_min) max(problem.u_max)];
        end
    end
    if ~isfield(opt,'y_lim') 
        opt.y_lim = [min(problem.y_min) max(problem.y_max)];
    end
    if ~isfield(opt,'yaxis_lim_enhance_rel')
        opt.yaxis_lim_enhance_rel = 1e-1; 
    end
    if ~isfield(opt,'linewidth_u')
        opt.linewidth_u = 1; 
    end
    if ~isfield(opt,'linewidth_y')
        opt.linewidth_y = 1; 
    end
    if ~isfield(opt,'linewidth_x')
        opt.linewidth_x = 1; 
    end
    if ~isfield(opt,'normalized_u') 
        opt.normalized_u = 0;
    end
    if ~isfield(opt,'normalized_y') 
        opt.normalized_y = 0;
    end
    if ~isfield(opt,'grid_on')
        opt.grid_on = 0;
    end
    
        
    %% initialization
    if problem.nU == 0
        opt.uy_same_axis = 1;
    end    
    if isempty(opt.idx_samples)
        opt.idx_samples = 1:length(dataset);
    end
    nS = length(opt.idx_samples);
    if opt.num_cols <= 0 || opt.num_rows <= 0
        opt.num_rows = ceil(sqrt(nS));
        opt.num_cols = ceil(nS/opt.num_rows);
    end
    if opt.yaxis_lim_enhance_rel > 0
        if problem.nU > 0
            opt.u_lim = opt.u_lim + [-1 1]*opt.yaxis_lim_enhance_rel*(opt.u_lim(2)-opt.u_lim(1));
        end
        opt.y_lim = opt.y_lim + [-1 1]*opt.yaxis_lim_enhance_rel*(opt.y_lim(2)-opt.y_lim(1));
    end
    
    for iS = 1:nS
        if isfield(dataset{iS},'uu')
            if isa(dataset{iS}.uu,'function_handle')
                if length(dataset{iS}.tt) == 2
                    dataset{iS}.tt = linspace(dataset{iS}.tt(1),dataset{iS}.tt(2),1e3);
                    dataset{iS}.uu = dataset{iS}.uu(dataset{iS}.tt);
                end
            end
        end
    end
     
    if opt.automargin && opt.plot_legend
       error('automargin and plot_legend incompatible') 
    end
    if opt.one_sample_per_time
        %% one_sample_per_time
        figure()
        for iS = 1:nS
            iSloc = opt.idx_samples(iS);
            subplot(2,1,1)
            
            if  opt.normalized_u
                plot(dataset{iSloc}.tt, (dataset{iSloc}.uu - problem.u_min)./(problem.u_max - problem.u_min), ...
                    'linewidth',opt.linewidth_u)
                axis([ min(dataset{iSloc}.tt) max(dataset{iSloc}.tt) 0 1 ])
            else
                plot(dataset{iSloc}.tt, dataset{iSloc}.uu, ...
                    'linewidth',opt.linewidth_u)
                axis([ min(dataset{iSloc}.tt) max(dataset{iSloc}.tt) min(problem.u_min) max(problem.u_max) ])
            end
                        
            title(sprintf('sample %d of %d',iS,nS))

            subplot(2,1,2)
            hold off
            if isfield(dataset{iSloc},'yy_ex')
                if  opt.normalized_y
                    plot(dataset{iSloc}.tt, (dataset{iSloc}.yy_ex - problem.y_min)./(problem.y_max - problem.y_min), ...
                    '-k', 'linewidth',opt.linewidth_y)
                    hold on
                    axis([ min(dataset{iSloc}.tt) max(dataset{iSloc}.tt) 0 1 ])
                else
                    plot(dataset{iSloc}.tt, dataset{iSloc}.yy_ex, ...
                    '-k', 'linewidth',opt.linewidth_y)
                    hold on
                    axis([ min(dataset{iSloc}.tt) max(dataset{iSloc}.tt) min(problem.y_min) max(problem.y_max) ])
                end                
            end
            if  opt.normalized_y
                plot(dataset{iSloc}.tt, (dataset{iSloc}.yy - problem.y_min)./(problem.y_max - problem.y_min), ...
                    'linewidth',opt.linewidth_y)
                axis([ min(dataset{iSloc}.tt) max(dataset{iSloc}.tt) 0 1 ])
            else
                plot(dataset{iSloc}.tt, dataset{iSloc}.yy, ...
                    'linewidth',opt.linewidth_y)
                axis([ min(dataset{iSloc}.tt) max(dataset{iSloc}.tt) min(problem.y_min) max(problem.y_max) ])
            end
            pause();
        end
    else
        %% single windows
        
        if opt.automargin            
            if isempty(opt.figure_size)
                figure()
            else
                figure('units','pixel','position',opt.figure_size);
            end
        else
            opt.size_opt.plot_legend = opt.plot_legend;
            win_sizes = plotting_get_axes(opt.num_rows,opt.num_cols,opt.size_opt);
            figure('units','pixel','position',[100 100 win_sizes.width win_sizes.heigth]);
        end
                
        for iS = 1:nS

            iSloc = opt.idx_samples(iS);
            
            if ~isfield(dataset{iSloc},'tt_y')
                dataset{iSloc}.tt_y = dataset{iSloc}.tt;
            end

            if opt.automargin
                subplot(opt.num_rows,opt.num_cols,iS)
            else
                subplot('Position',win_sizes.axs{iS}.coord_norm)
            end
            
            if problem.nU > 0
                if ~opt.uy_same_axis
                    yyaxis left
                end
                if  opt.normalized_u 
                    plot(dataset{iSloc}.tt,(dataset{iSloc}.uu - problem.u_min)./(problem.u_max - problem.u_min), ...
                        '-','linewidth',opt.linewidth_u); hold on;
                    axis([0 max(dataset{iSloc}.tt) 0 1]) 
                else
                    plot(dataset{iSloc}.tt,dataset{iSloc}.uu, ...
                        '-','linewidth',opt.linewidth_u); hold on;
                    axis([0 max(dataset{iSloc}.tt) opt.u_lim]) 
                end
            end

            if ~opt.uy_same_axis
                yyaxis right
            end
            if isfield(dataset{iSloc},'yy_ex')
                if  opt.normalized_y
                    plot(dataset{iSloc}.tt,(dataset{iSloc}.yy_ex - problem.y_min)./(problem.y_max - problem.y_min), ...
                        'k-','linewidth',opt.linewidth_y); hold on;
                else
                    plot(dataset{iSloc}.tt,dataset{iSloc}.yy_ex, ...
                        'k-','linewidth',opt.linewidth_y); hold on;
                end
            end
            if isfield(dataset{iSloc},'yy')
                if  opt.normalized_y
                    plot(dataset{iSloc}.tt_y,(dataset{iSloc}.yy - problem.y_min)./(problem.y_max - problem.y_min), ...
                        '-','linewidth',opt.linewidth_y); hold on;
                else
                    plot(dataset{iSloc}.tt_y,dataset{iSloc}.yy, ...
                        '-','linewidth',opt.linewidth_y); hold on;
                end
            end
            if  opt.normalized_y
                axis([0 max(dataset{iSloc}.tt) 0 1]) 
            else
                axis([0 max(dataset{iSloc}.tt) opt.y_lim]) 
            end
            if opt.grid_on
                grid on
            end
            
            if opt.plot_x
                if isequal(opt.plot_x_idxs, -1)
                    opt.plot_x_idxs = 1:size(dataset{iSloc}.xx,1);
                end
                    
                plot(dataset{iSloc}.tt, dataset{iSloc}.xx(opt.plot_x_idxs,:), ...
                    'k--','linewidth',opt.linewidth_x); hold on;
            end
            
%             Tmax = max(set{iSloc}.tt);
% 
%             plot(set{iSloc}.tt,(set{iSloc}.uu(1,:)-opts.lim_u(1,1))/(opts.lim_u(1,2)-opts.lim_u(1,1)),'linewidth',linewidthplots,'Color',col_u1); hold on;
%             if opts.sarcomeres
%                 yyaxis right
%                 plot(set{iSloc}.tt,set{iSloc}.uu(2,:),'linewidth',linewidthplots,'Color',col_u2); hold on;
%                 axis([0 Tmax opts.u2Limits])
%                 ax = gca;
%                 ax.YColor = [0 0 0];
%                 yyaxis left
%             else
%                 plot(set{iSloc}.tt,(set{iSloc}.uu(2,:)-opts.lim_u(2,1))/(opts.lim_u(2,2)-opts.lim_u(2,1)),'linewidth',linewidthplots,'Color',col_u2); hold on;
%             end
% 
%             plot(set{iSloc}.tt,(set{iSloc}.yy-opts.lim_y(1))/(opts.lim_y(2)-opts.lim_y(1)),'linewidth',linewidthplots,'Color',col_y); hold on;
% 
%             axis([0 Tmax opts.yLimits])   

            if isfield(dataset{iSloc},'label')
                title(dataset{iSloc}.label)
            end

            if opt.show_x_ticks == 0
                set(gca,'XTickLabel',[]);
            elseif opt.show_x_ticks == -1
                if win_sizes.axs{iS}.row < opt.num_rows
                    set(gca,'XTickLabel',[]);
                else
                    xlabel('t')
                end
            end
            if opt.show_y_ticks == 0
                if problem.nU > 0
                    if ~opt.uy_same_axis, yyaxis left, end
                    set(gca,'YTickLabel',[]);
                    ax = gca;
                    ax.YColor = [0 0 0];
                    if ~opt.uy_same_axis
                        yyaxis right
                        set(gca,'YTickLabel',[]);
                        ax = gca;
                        ax.YColor = [0 0 0];
                    end
                else
                    set(gca,'YTickLabel',[]);
                end
            elseif opt.show_y_ticks == -1
                if problem.nU > 0
                    if ~opt.uy_same_axis
                        if win_sizes.axs{iS}.col < opt.num_cols
                            yyaxis right
                            set(gca,'YTickLabel',[]);
                        else
                            yyaxis right
                            ylabel(opt.y_label)
                        end
                    end
                    if win_sizes.axs{iS}.col > 1
                        if ~opt.uy_same_axis, yyaxis left, end
                        set(gca,'YTickLabel',[]);
                    else
                        if ~opt.uy_same_axis
                            yyaxis left
                            ylabel('u')
                        end
                    end
                else
                    if win_sizes.axs{iS}.col > 1
                        set(gca,'YTickLabel',[]);
                    else
                        ylabel(opt.y_label)
                    end
                end
            end
            
            pause(1e-16)

        end
                    
        if opt.plot_legend
            myleg = legend(problem_getvariablename_list(problem),'Orientation','horizontal');
            set(myleg,'position', win_sizes.get_legend_coord_norm(myleg.Position), 'units', 'normalized');           
        end

        
    end
    if opt.do_save
        print(opt.filename,'-depsc','-painters');
    end
    
    if opt.do_save && opt.do_save_external_legend
        iLine = 1;
        for i = 1:problem.nU
            lines{iLine}.width = opt.linewidth_u;
            lines{iLine}.name = problem_getvariablename(problem,'u',i);
%             if isfield(problem,'u_names')
%                 lines{iLine}.name = problem.u_names{iLine};
%             else
%                 if problem.nU == 1
%                     lines{iLine}.name = 'u';
%                 else
%                     lines{iLine}.name = sprintf('u_%d',i);
%                 end
%             end
            iLine = iLine + 1;
        end
        for i = 1:problem.nY
            lines{iLine}.width = opt.linewidth_y;
            lines{iLine}.name = problem_getvariablename(problem,'y',i);
%             if isfield(problem,'y_names')
%                 lines{iLine}.name = problem.y_names{iLine};
%             else
%                 if problem.nY == 1
%                     lines{iLine}.name = 'y';
%                 else
%                     lines{iLine}.name = sprintf('y_%d',i);
%                 end
%             end
            iLine = iLine + 1;
        end
        dataset_plot_legend(lines,[opt.filename '_legend'])
    end
    fprintf(' done!\n')
end

