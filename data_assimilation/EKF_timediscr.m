function out = EKF_timediscr(f,g,dfdx,dgdx,x0,P0,Q,R,yy_ex,uu,opt)

% Extended Kalman Filter - discrete time case
%
% x^{k+1} = f(x^k,u^k,u^{k+1}) + w^k      w^k: N(0,Q),   Q:  [X^2]
% y^k = g(x^k)+ v^k                       v^k: N(0,R),   R:  [Y^2]
%
% E[x^1] = x0, Var[x^{1}] = P0                           P0: [X^2]

    opt.dummy = 0;
    if ~isfield(opt,'plot_freq') % -1 = never; 0 = last time interval; n = every n step
        opt.plot_freq = -1;
    end
    if ~isfield(opt,'idx_param')
        opt.idx_param = [];
    end
    if ~isfield(opt,'export_x')
        opt.export_x = 1;
    end
    if ~isfield(opt,'export_y')
        opt.export_y = 1;
    end
    if ~isfield(opt,'export_Pdiag')
        opt.export_Pdiag = 1;
    end
    if ~isfield(opt,'do_tests')
        opt.do_tests = 0;
    end
    
    nT = size(yy_ex,2);
    nX = size(x0,1);
    nY = size(yy_ex,1);
    
    x = x0;
    P = P0;
    
    if isempty(dfdx)
        dfdx = @(x,u1,u2) numerical_differentiation(@(z) f(z,u1,u2),x);
        fprintf('warning: using numerical differentiation for dfdx!\n')
    end
    if isempty(dgdx)
        dgdx = @(x) numerical_differentiation(g,x);
        fprintf('warning: using numerical differentiation for dgdx!\n')
    end
    
    if opt.do_tests
        perturbation = @(v) v + min(1e-8, 1e-2*abs(v)).*rand(size(v));
        x_base = perturbation(x0);
        u1_base = perturbation(uu(:,1));
        u2_base = perturbation(uu(:,2));
        
        f_test = @(x) f(x,u1_base,u2_base);
        
        dfdx_teo = dfdx(x_base,u1_base,u2_base);
        dfdx_num = numerical_differentiation(f_test,x_base);
        dgdx_teo = dgdx(x_base);
        dgdx_num = numerical_differentiation(g,x_base);
        
        fprintf('numerical differentiation test:\n')
        fprintf('dfdx %1.2e\n',norm(dfdx_num(:)-dfdx_teo(:))/norm(dfdx_teo(:)))
        fprintf('dgdx %1.2e\n',norm(dgdx_num(:)-dgdx_teo(:))/norm(dgdx_teo(:)))
        
%         numerical_diff_testconvergence(f_test,x0,struct('jac_exact',dfdx_teo))
%         numerical_diff_testconvergence(g,x0,struct('jac_exact',dgdx_teo))
%         opt.idx_out = [];
%         opt.idx_in = [];
    end
        
    if opt.export_x
        xx = nan(nX,nT);
        xx_pred = nan(nX,nT);
        xx(:,1) = x;
    end
    if opt.export_y
        yy = nan(nY,nT);
        yy_pred = nan(nY,nT);
        yy(:,1) = g(x);
    end
    if opt.export_Pdiag
        PP = nan(nX,nT);
        PP(:,1) = diag(P);
    end
    
    if opt.plot_freq >= 0
%         tt = dt*(0:nT-1);
        tt = 1:nT;
        FigMain = figure('units','normalized','outerposition',[0 0 1 1]);
        cmap = get(0, 'DefaultAxesColorOrder');
        nrows = 2;
        
        idx_param = opt.idx_param;
        idx_state = setdiff(1:nX,idx_param);
        
        if ~isempty(idx_param)
            nrows = nrows+1;
        end
    end
    
    for iT = 2:nT
    
        F = dfdx(x,uu(:,iT-1),uu(:,iT));
        
        % Prediction
        x_pred = f(x,uu(:,iT-1),uu(:,iT)); % NB: explicit in u
        P_pred = F*P*F' + Q;
        y_pred = g(x_pred);
        
        % Correction
        e = yy_ex(:,iT) - y_pred;
        H = dgdx(x_pred);
        C = H*P_pred*H' + R;
        K = P_pred*H'/C;
        x = x_pred + K*e;
        P = P_pred - K*H*P_pred;
        
        if opt.export_x
            xx(:,iT) = x;
            xx_pred(:,iT) = x_pred;
        end
        if opt.export_y
            yy(:,iT) = g(x);
            yy_pred(:,iT) = y_pred;
        end
        if opt.export_Pdiag
            PP(:,iT) = diag(P);
        end
        
        if isfield(opt, 'postprocessing')
            opt.postprocessing(iT,x);
        end
        
        if opt.plot_freq >= 0 
            if mod(iT,opt.plot_freq) == 0 || iT == nT
                figure(FigMain)
                if opt.export_y 
                    subplot(nrows,1,1)
                    hold off
                    plot(tt(1:iT),yy_ex(:,1:iT),'--'); hold on
                    plot(tt(1:iT),yy_pred(:,1:iT)); hold on
                end             
                if opt.export_x   
         
                    subplot(nrows,1,2)
                    hold off
                    plot(tt(1:iT),xx(idx_state,1:iT),'linewidth',1.5); hold on
                    
                    if ~isempty(idx_param)
                        subplot(nrows,1,3)
                        hold off
                        plot(tt(1:iT),xx(idx_param,1:iT),'linewidth',1.5); hold on
                    end
                    
                    if opt.export_Pdiag
                        band = 3*sqrt(PP(:,1:iT));
                        xx_up = xx(:,1:iT) + band;
                        xx_dw = xx(:,1:iT) - band; 
                        
                        subplot(nrows,1,2)
                        for iX=1:length(idx_state)
                            col = cmap(1 + mod(iX-1,size(cmap,1)),:);
                            patch([tt(1:iT) tt(iT:-1:1)],[xx_dw(idx_state(iX),:) xx_up(idx_state(iX),end:-1:1)],col,'EdgeColor','none'); hold on
                            alpha(.3) 
                        end
                        
                        if ~isempty(idx_param)
                            subplot(nrows,1,3)
                            for iX=1:length(idx_param)
                                col = cmap(1 + mod(iX-1,size(cmap,1)),:);
                                patch([tt(1:iT) tt(iT:-1:1)],[xx_dw(idx_param(iX),:) xx_up(idx_param(iX),end:-1:1)],col,'EdgeColor','none'); hold on
                                alpha(.3) 
                            end
                        end
                    end
                    
%                     subplot(nrows,1,2)
%                     ax = gca;
%                     ax.ColorOrderIndex = 1;
%                     plot(tt(1:iT),xx(idx_state,1:iT),'linewidth',1.5); hold off
%                     
%                     if ~isempty(idx_param)
%                         subplot(nrows,1,3)
%                         ax = gca;
%                         ax.ColorOrderIndex = 1;
%                         plot(tt(1:iT),xx(idx_param,1:iT),'linewidth',1.5); hold off
%                     end
                    
                    
                    subplot(nrows,1,2)
                    for iX=1:length(idx_state)
                        col = cmap(1 + mod(iX-1,size(cmap,1)),:);
                        plot(tt(1:iT),xx(idx_state(iX),1:iT),'Color',col,'linewidth',1.5); hold on
                    end
                    hold off

                    if ~isempty(idx_param)
                        subplot(nrows,1,3)
                        for iX=1:length(idx_param)
                            col = cmap(1 + mod(iX-1,size(cmap,1)),:);
                            plot(tt(1:iT),xx(idx_param(iX),1:iT),'Color',col,'linewidth',1.5); hold on
                        end
                        hold off
                    end
                    
                end
                pause(1e-16)
            end
        end
        
    end
    
    out.x = x;
    out.P = P;
%     out.idx_param = idx_param;
    if opt.export_x
        out.xx      = xx;
        out.xx_pred = xx_pred;
    end
    if opt.export_y
        out.yy      = yy;
        out.yy_pred = yy_pred;
    end
    if opt.export_Pdiag
        out.PP = PP;
    end

end