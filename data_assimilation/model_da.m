function out_da = model_da(metamod,test,opt)

    %% Options
    opt.dummy = 0;
    if ~isfield(opt,'obs_err') % [Y*sqrt(T)]
        opt.obs_err = 1e-1;
    end
    if ~isfield(opt,'mod_err') % [X/sqrt(T)]
        opt.mod_err = 0;
    end
    if ~isfield(opt,'do_plot')
        opt.do_plot = 0;
    end
    if ~isfield(opt,'numerical_diff')
        opt.numerical_diff = 0;
    end
    
    %% Initialization
    N = metamod.nX;
    if isfield(metamod, 'nA')
        nA = metamod.nA;
    else
        nA = 0;
    end
    nY = metamod.problem.nY;
    nU = metamod.problem.nU;
    
    if ~isfield(test,'tt_y')
        test.tt_y = test.tt;
    end
    
    metamod = adapt_dimension_struct_field(metamod,'x0_min',N);
    metamod = adapt_dimension_struct_field(metamod,'x0_max',N);
    metamod = adapt_dimension_struct_field(metamod,'alpha_min',nA);
    metamod = adapt_dimension_struct_field(metamod,'alpha_max',nA);
    opt = adapt_dimension_struct_field(opt,'obs_err',nY);
    opt = adapt_dimension_struct_field(opt,'mod_err',N);

    % if the model started in the past, use the generic apriori for x0
    if test.tt_y(1) == test.tt(1)
        [x0_mean,x0_wide] = metamodel_get_stat_x0(metamod);
    else
        [x0_mean,x0_wide] = metamodel_get_stat_x(metamod);
    end
    x0 = x0_mean;
    P0 = diag(x0_wide.^2);
    if nA > 0
        [a0_mean,a0_wide] = metamodel_get_stat_alpha(metamod);
        x0 = [x0;a0_mean];
        P0(N+1:N+nA,N+1:N+nA) = diag(a0_wide.^2);
    end
    Q = diag([opt.mod_err; zeros(nA,1)].^2);
    R = diag(opt.obs_err.^2);
    dt = metamod.dt;

    % advance functions
    if nA > 0
        switch metamod.advance_type
            case 'nonlinear_explicit'
                f = @(x,u) [metamod.f_alpha(x(1:N),u,x(N+1:end)); zeros(nA,1)];
                dfdx = @(x,u) [metamod.dfdx(x(1:N),u,x(N+1:end)), metamod.dfda(x(1:N),u,x(N+1:end)); zeros(nA,N+nA)];
                f_discr = @(x,u1,u2) x+dt*f(x,u1);
                dfdx_discr = @(x,u1,u2) eye(size(x0,1)) + dt*dfdx(x,u1);
%             case 'linear_advance_timeexplicit_uaffine'   
%                 f_discr = @(x,u1,u2) [xnew_lateua(x(1:N),u1,x(N+1:end),dt); zeros(nA,1)];
%                 dfdx_discr = @(x,u1,u2) [dxnewdx_lateua(x(1:N),u1,x(N+1:end),dt), dxnewda_lateua(x(1:N),u1,x(N+1:end),dt); zeros(nA,N+nA)];
            case 'linear_CN_uaffine'
                f_discr = @(x,u1,u2) [xnew_CN(x(1:N),u1,x(N+1:end),dt); x(N+1:end)];
                dxnewdx = @(a) (eye(N)/dt - metamod.A_alpha(a)/2)\(eye(N)/dt + metamod.A_alpha(a)/2);
                dfdx_discr = @(x,u1,u2) [dxnewdx(x(N+1:end)), dxnewda_CN(x(1:N),u1,x(N+1:end),dt); zeros(nA,N), eye(nA)];            
            otherwise
                error('not (yet) implemented')
        end
    else
        switch metamod.advance_type
            case 'nonlinear_explicit'
                dfdx = @(x,u) metamod.dfdx(x,u,[]);
                f_discr = @(x,u1,u2) x+dt*metamod.f(x,u1);
                dfdx_discr = @(x,u1,u2) eye(size(x0,1)) + dt*dfdx(x,u1);     
            otherwise
                error('not (yet) implemented')
        end
    end

    % output functions
    if isequal(metamod.output_type,'insidestate')
        g = @(x) x(1:nY);
        dgdx = @(x) [eye(nY), zeros(nY,N+nA-nY)];
    elseif isequal(metamod.output_type,'nonlinear')
        g = @(x) metamod.g(x(1:N));
        dgdx = @(x) [metamod.dgdx(x(1:N)), zeros(nY,nA)];
    elseif isequal(metamod.output_type,'linear')
        g = @(x) metamod.G*x(1:N)+metamod.g0;
        dgdx = @(x) [metamod.G, zeros(nY,nA)];
    end
    
    % fixing noise homogeneity
    Q_discr = dt*Q;
    R_discr = R/dt;
    
    if nU == 0
        test.uu = zeros(1,size(test.yy,2));
    end
    
    %% Execution
    tt = test.tt_y(1):dt:test.tt_y(end);
    % TODO: instead of interpolating one should account for observation at isolate timesteps
    uu = interp_time_series(test.tt,test.uu,tt);
    yy = interp_time_series(test.tt_y,test.yy,tt);
%     out_da = EKF_timect_FE(f,g,dfdx,dgdx,x0,P0,Q,R,yy,uu,dt);
%     if opt.do_plot
%         da_opt.plot_freq = 0;
%     end
    da_opt.dummy = 0;
    da_opt.do_tests = 0;
    da_opt.idx_param = N+1:N+nA;
    if opt.numerical_diff
        dfdx_discr = [];
        dgdx = [];
    end
    out_da = EKF_timediscr(f_discr,g,dfdx_discr,dgdx,x0,P0,Q_discr,R_discr,yy,uu,da_opt);
    out_da.tt = tt;
    
    %% plot    
    if opt.do_plot
        da_plot_output(out_da,test,metamod);
%         linewdt = 1.2;
%         figure();
%         nrows = 1;
%         ncols = 1;
%         if metamod.nA > 0
%             nrows = 2; 
%         end
%         subplot(nrows,ncols,1)
%         plot(test.tt,test.yy,'-','linewidth',linewdt); hold on; ax = gca; ax.ColorOrderIndex = 1;
% %         plot(test.tt,out_HF_noised.yy,'-','linewidth',1); hold on; ax = gca; ax.ColorOrderIndex = 1;
%         if isequal(metamod.output_type,'insidestate')
%             plot(out_da.tt,out_da.xx(1:metamod.problem.nY,:),':','linewidth',linewdt); hold on; ax = gca; ax.ColorOrderIndex = 1;
%             plot(out_da.tt(2:end),out_da.xx_pred(1:metamod.problem.nY,2:end),'--','linewidth',linewdt); hold on; ax = gca; ax.ColorOrderIndex = 1;
%         else
%             plot(out_da.tt(2:end),out_da.yy_pred(:,2:end),'--','linewidth',linewdt); hold on; ax = gca; ax.ColorOrderIndex = 1;
%         end
%         % axis([0 Tmax 0 3])
%         if metamod.nA > 0
%             subplot(nrows,ncols,2)
%     %         plot(out_HF.tt,out_da.xx(metamod.problem.nY+1:end,:),':','linewidth',linewdt); hold on; ax = gca; ax.ColorOrderIndex = 1;
%     %         plot(out_HF.tt(2:end),out_da.xx_pred(metamod.problem.nY+1:end,2:end),'--','linewidth',linewdt); hold on; ax = gca; ax.ColorOrderIndex = 1;
%             plot(out_da.tt,out_da.xx(metamod.nX+1:end,:),':','linewidth',linewdt); hold on; ax = gca; ax.ColorOrderIndex = 1;
%             plot(out_da.tt(2:end),out_da.xx_pred(metamod.nX+1:end,2:end),'--','linewidth',linewdt); hold on; ax = gca; ax.ColorOrderIndex = 1;
% %             plot(out_da.tt([1 end]),mod.alpha*[1 1])
%         end
%         % axis([0 Tmax a_min a_max])
%         pause(1e-16)
    end
    
    %% Functions for particual advance types:
    
    %% --- Crank-Nicolson
    
    function xnew = xnew_CN(x,u,a,dt)
%         if isfield(metamod,'b0_alpha')
%             b  = metamod.b0_alpha(a);
%         else
%             b  = metamod.b0;
%         end
%         for iU = 1:nU
%             if isfield(metamod,'b0_alpha')
%                 b  = b  + u(iU)*metamod.b_alpha{iU}(a);
%             else
%                 b  = b  + u(iU)*metamod.b{iU};
%             end
%         end
%         if isfield(metamod,'A_alpha')
%             A  = metamod.A_alpha(a);
%         else
%             A  = metamod.A;
%         end    
        b  = metamod.b0_alpha(a);
        for iU = 1:nU
            b  = b  + u(iU)*metamod.b_alpha{iU}(a);
        end
        A  = metamod.A_alpha(a);
        xnew = (eye(N)/dt - .5*A) \ ((eye(N)/dt + .5*A)*x + b); 
    end

    function dxnewda = dxnewda_CN(x,u,a,dt)
        x2 = xnew_CN(x,u,a,dt);
        dxnewda = (eye(N)/dt - metamod.A_alpha(a)/2)\((dhdw(x2,a,u)+dhdw(x,a,u))/2);
    end
    
    function D = dhdw(x,a,u)
        D  = metamod.db0_dalpha(a);
        for iU = 1:nU
            D  = D  + u(iU)*metamod.db_dalpha{iU}(a);
        end
        dAdalpha = metamod.dA_dalpha(a);
        for iW = 1:nA
            D(:,iW) = D(:,iW) + dAdalpha(:,:,iW)*x;
        end
    end

    %% --- linear_advance_timeexplicit_uaffine
    
%     function out = xnew_lateua(x,u,a,dt,outtype)
%         Kf = metamod.Kf0;
%         Kv = metamod.Kv0;
%         Bf = metamod.Bf0;
%         Bv = metamod.Bv0;
%         b  = metamod.b0;
%         for iU = 1:nU
%             Kf = Kf + u(iU)*metamod.Kf{iU};
%             Kv = Kv + u(iU)*metamod.Kv{iU};
%             Bf = Bf + u(iU)*metamod.Bf{iU};
%             Bv = Bv + u(iU)*metamod.Bv{iU};
%             b  = b  + u(iU)*metamod.b{iU};
%         end
%         switch outtype
%             case 0
%                 out = (Kf + Kv/dt) \ ((Bf + Bv/dt)*x + b); 
%             case 1
%                 out = (Kf + Kv/dt) \ (Bf + Bv/dt); 
%             case 2
%         end
%     end

end