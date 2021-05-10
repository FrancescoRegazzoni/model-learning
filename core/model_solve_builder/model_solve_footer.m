        timeElapsed = toc(timeInit);

        %% postprocessing
        output.tt = tt;
        if nU > 0
            output.uu = uu;
        end
        output.tt_y = tt_y;
        output.yy = yy;
        output.time = timeElapsed;
        output.time_norm = timeElapsed/(tt(end)-tt(1));
        if opt.save_x
            output.xx = xx;
        end
        if opt.error_compute
            output.err_L2 = get_norm_L2_time(tt_y,yy-yy_ex);
    %         output.err_L2 = sqrt(mean((yy(:) - yy_ex(:)).^2));
        end
        if isfield(test,'yy') && opt.save_y_ex
            output.yy_ex = yy_ex;
        end
        if isfield(test,'label')
            output.label = test.label;
        end

        if opt.verbose
            fprintf('   model solved --- time elapsed: %1.2e s', timeElapsed)
            if opt.error_compute
                fprintf(' --- error %1.2e',output.err_L2);
            end
            fprintf('\n');
        end

    end
                    
    if opt.do_plot
        nrows = 1;
        if model.problem.nU > 0
            nrows = nrows+1;
        end
        if opt.do_plot_x 
            nrows = nrows+1;
        end
        
        i_row = 1;
        
        if model.problem.nU > 0
            subplot(nrows,1,i_row)
            plot(output.tt,output.uu,'-','linewidth',1.2)
            axis([0 output.tt(end) min(model.problem.u_min) max(model.problem.u_max)])
            ylabel('u')
            i_row = i_row + 1;
        end
        
        subplot(nrows,1,i_row)
        plot(output.tt_y,output.yy,'-','linewidth',1.2) 
        axis([0 output.tt(end) min(model.problem.y_min) max(model.problem.y_max)])
        if isfield(test,'yy')
            hold on
            sty = '--'; lnwdt = 1.2;
            if nY > 1
                % if the number of outputs ig greater than 1, I use the same colour
                ax = gca;
                ax.ColorOrderIndex = 1;
                sty = '--'; lnwdt = 1.5;
            end
            plot(tt_y,yy_ex,sty,'linewidth',lnwdt)
            if nY == 1
                legend('test model','HF model')
            end
            hold off
        end
        ylabel('y')
        i_row = i_row + 1;
        
        if opt.do_plot_x && ~model.blackbox
            subplot(nrows,1,i_row)
            if strcmp(model.advance_type,'ANN')
                plot(tt,xx(nY+1:end,:),'-','linewidth',1)
            else
                plot(tt,xx,'-','linewidth',1)
            end
            ylabel('x')
        end
        
        pause(1e-16)
    end
end