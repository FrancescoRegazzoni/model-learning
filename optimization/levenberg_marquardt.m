function x0 = levenberg_marquardt(func,x0,opt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Consider the following problem
%       min f(x) = 1/2 \sum_{i=1}^n r_i(x), where x \in \mathbb{R}^N
%
% x0            Initial guess
% FuncGrad      Function with the following interface:
%               INPUTS:   x (N,1)
%               OUTPUTS:  f (1,1)
%                         df/dx (N,1)
% LSFuncGrad    Function with the following interface:
%               INPUTS:   x (N,1)
%               OUTPUTS:  R (n,1)
%                         dR/dx (n,N)
% opt           options structures. Fields:
%               - nMaxIter
%               - alphaFix
%               - SA
%               - SA_stepsize
%               - SA_T
%               - LSopt.alphaFix         
%               - LSopt.algorithm         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% options
    opt.dummy = 0;
    opt.LSopt.dummy = 0;    
    if ~isfield(opt,'nMaxIter')
        opt.nMaxIter = 1e6;
    end
    if ~isfield(opt,'do_plot')
        opt.do_plot = 0;
    end
%     if ~isfield(opt,'verbose')
%         opt.verbose = 1;
%     end
    if ~isfield(opt,'tol_increment')
        opt.tol_increment = 1e-12;
    end
    if ~isfield(opt,'dphiZeroThreshold')
        opt.dphiZeroThreshold = -1;
    end
    if ~isfield(opt,'SA')
        opt.SA = 0;
    end
    if ~isfield(opt.LSopt,'algorithm')
        opt.LSopt.algorithm = 2;
    end
    if ~isfield(opt.LSopt,'alphaFix')
        opt.LSopt.alphaFix = .1;
    end
    
    %% input functions
    if ~isfield(func,'residuals') && ~isfield(func,'residuals_jacobian')
        error('Levenberg-Marquardt: function to evaluate residuals must be present')
    end
    numerical_jac = ~isfield(func,'residuals_jacobian') && ~isfield(func,'jacobian');
    if numerical_jac
        func.residuals_jacobian = @residuals_jacobian;
    elseif ~isfield(func,'residuals_jacobian') && isfield(func,'jacobian')        
        func.residuals_jacobian = @residuals_jacobian_from_res_and_jac;
    end
    if ~isfield(func,'cost_gradient')
        if numerical_jac
            if ~isfield(func,'cost')
                func.cost = @(x) .5 * sum(func.residuals(x).^2);
            end            
            func.cost_gradient = @cost_gradient_numerical;
        else
            func.cost_gradient = @cost_gradient_fromjacobian;
        end
    end
    
    if ~isfield(func,'step_finalization')
        func.step_finalization = @(x) x;
    end
       
    function [cost,grad] = cost_gradient_numerical(x)
        cost = func.cost(x);
        if nargout == 2
            grad = numerical_diff(func.cost,x,struct('F0',cost))';
        end
    end     
    
    function [cost,grad] = cost_gradient_fromjacobian(x)        
        if nargout == 1
            res = func.residuals_jacobian(x);
        else
            [res,jac] = func.residuals_jacobian(x);
        end
        cost = .5 * sum(res.^2);
        if nargout == 2
            grad = jac'*res;
        end
    end     
    
    function [res,jac] = residuals_jacobian(x)
        res = func.residuals(x);
        if nargout == 2
            jac = numerical_diff(func.residuals,x,struct('F0',res));
        end
    end

    function [res,jac] = residuals_jacobian_from_res_and_jac(x)
        res = func.residuals(x);
        if nargout == 2
            jac = func.jacobian(x);
        end
    end

    %% execution    
    it = 0;
    n = length(x0);
    x0 = func.step_finalization(x0);
    mu1 = 1;
        
    E0 = func.cost_gradient(x0);
    mu_ref = mu1*E0;
        
    opt.LSopt.dphiZeroThreshold = opt.dphiZeroThreshold;
    
    if opt.do_plot
        fig_out = figure();
        E_hist = E0;
        E_diff_hist = [];
        x_diff_hist = [];
        x_old = x0;
    end
    
    if opt.tol_increment > 0 || opt.do_plot
        x_old = x0;
    end
    
    goon = 1;
    
    while (it < opt.nMaxIter && goon)
        TimeIterInit = tic();
        
        it = it+1;
        fprintf('LM iteration %d',it)
        
        TimeJacobianInit = tic();
        [F,DF] = func.residuals_jacobian(x0);
        TimeJacobian = toc(TimeJacobianInit);
        
        TimeLMInit = tic();
%         mu = min(mu1,norm(DF'*F));
%         mu = min(F'*F,norm(DF'*F))
%         mu = min(mu_ref,norm(DF'*F));

        mu_GS = norm(DF'*F); % Grippo-Sciandrone
        if mu_ref < mu_GS
            mu = mu_ref;
            fprintf(' --- constant mu')
        else
            mu = mu_GS;
            fprintf(' --- GS mu')
        end
        fprintf(' ( mu = %1.2e )\n',mu)
            
        dx = (DF'*DF + mu*eye(n))\(DF'*F);
        TimeLM = toc(TimeLMInit);

        TimeLSInit = tic();
        
        phiZero = .5*(F'*F);
        dphiZero = (DF'*F)'*(-dx);
    
        if opt.LSopt.algorithm == 0
            x0 = x0 - opt.LSopt.alphaFix*dx;
        else
            [x0,alphaStar] = linearsearch_vectorial(x0,-dx,func.cost_gradient,opt.LSopt,phiZero,dphiZero);
            opt.LSopt.alpha0 = alphaStar;
        end
        TimeLS = toc(TimeLSInit);
        
        if opt.tol_increment > 0 || opt.do_plot
            incr = norm(x0-x_old);
            x_old = x0;
            if incr < opt.tol_increment
                goon = 0;
            end
        end
        
        if abs(dphiZero/phiZero) < opt.dphiZeroThreshold
            fprintf('**************************************************************************************\n')
            fprintf('Activating Simulated Annealing...\n')
            fprintf('**************************************************************************************\n')
            opt.SA = 1;
            opt.SA_stepsize = @(k) 1e-4;
            opt.SA_T = @(k) 1;
        end
        
        TimeSAInit = tic();
        if isfield(opt,'SA')
            if opt.SA
                x0guess = x0 + randn(n,1) * opt.SA_stepsize(it);
                E0 = Func(x0);
                E0guess = func.cost_gradient(x0guess);
                accepted = 0;
                if rand(1) < exp(-(E0guess-E0)/E0/opt.SA_T(it))
                    accepted = 1;
                    x0 = x0guess;
                end
                fprintf('   SA: step = %1.2e  -  T = %1.2e  -  DeltaE = %1.2e  -  accepted = %d\n',opt.SA_stepsize(it),opt.SA_T(it),(E0guess-E0),accepted);
            end
        end
        TimeSA = toc(TimeSAInit);
                
        TimePlotInit = tic();
        x0 = func.step_finalization(x0);
        TimePlot = toc(TimePlotInit);
        
        TimeTot = toc(TimeIterInit);
        %PlotFnc(x0);
        fprintf('   Iter %d   Time: %1.2e   Time Jac: %1.2e   Time LM: %1.2e   Time LS: %1.2e   Time SA: %1.2e   Time plot: %1.2e \n', it, TimeTot, TimeJacobian, TimeLM, TimeLS, TimeSA, TimePlot)
                
        if opt.do_plot
            figure(fig_out);
            E_hist = [E_hist func.cost(x0)];
            E_diff_hist = [E_diff_hist abs(E_hist(end)-E_hist(end-1))];
            x_diff_hist = [x_diff_hist norm(x0-x_old)];
            x_old = x0;
            
            subplot(2,1,1)
            loglog(E_hist,'o-');
            subplot(2,1,2)
            loglog(E_diff_hist,'o-'); hold on
            loglog(x_diff_hist,'o-'); hold on
            hold off
            pause(1e-16)
        end
    end
    
    fprintf('Levenberg-Marquardt terminated after %d iterations\n',it)
end