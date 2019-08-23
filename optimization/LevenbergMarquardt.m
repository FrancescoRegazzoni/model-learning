function LevenbergMarquardt(x0,Func,FuncGrad,LSFuncGrad,PlotFnc,opt)
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
    it = 0;
    n = length(x0);
    x0 = PlotFnc(x0);
    mu1 = 1;
        
    E0 = Func(x0);
    mu_ref = mu1*E0;
    
    opt.dummy = 0;
    opt.LSopt.dummy = 0;
    
    if ~isfield(opt,'nMaxIter')
        opt.nMaxIter = 1e6;
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
    
    opt.LSopt.dphiZeroThreshold = opt.dphiZeroThreshold;    
    
    while (it < opt.nMaxIter)
        TimeIterInit = tic();
        
        it = it+1;
        fprintf('LM iteration %d',it)
        
        TimeJacobianInit = tic();
        [F,DF] = LSFuncGrad(x0);
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
            [x0,alphaStar] = LinearSearchVectorial(x0,-dx,Func,FuncGrad,opt.LSopt,phiZero,dphiZero);
            opt.LSopt.alpha0 = alphaStar;
        end
        TimeLS = toc(TimeLSInit);
        
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
                E0guess = Func(x0guess);
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
        x0 = PlotFnc(x0);
        TimePlot = toc(TimePlotInit);
        
        TimeTot = toc(TimeIterInit);
        %PlotFnc(x0);
        fprintf('   Iter %d   Time: %1.2e   Time Jac: %1.2e   Time LM: %1.2e   Time LS: %1.2e   Time SA: %1.2e   Time plot: %1.2e \n', it, TimeTot, TimeJacobian, TimeLM, TimeLS, TimeSA, TimePlot)
    end
end