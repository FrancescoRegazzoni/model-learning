function pen_handler = penalization_alpha_equilibria(optionsfile)

    Pen_s           = iniread(optionsfile,'Penalizations','pen_alpha_eq','s','0');
    u_star          = iniread(optionsfile,'Penalizations','pen_alpha_eq_ustar','d','0');
    
    pen_handler.num_pen = 1;
    pen_handler.names = {'ERRalphaEq'};
    pen_handler.coefficient_string = {Pen_s};
    pen_handler.dependence_f = 1;
    pen_handler.dependence_g = 0;
    pen_handler.dependence_a = 0;
    pen_handler.dependence_i = 0;
    
    AlphaGrid = [];
    N = [];
    N_alpha = [];
    nwF = [];
    nS = [];
    num_equilibria = [];
    useG = [];
       
    pen_handler.initialize = @initialize;
    function pen_handler = initialize(pen_handler,problem,N_,N_alpha_,useG_,nS_,nwF_,nwG,nAL,nIC,misc)
        if pen_handler.pen_active
            N = N_;
            N_alpha = N_alpha_;
            nwF = nwF_;
            nS = nS_;
            useG = useG_;
            normaliz = misc.normaliz;
            if mod(N_alpha, N) ~= 0
                error('penalization of alpha_eq works only if N_alpha is a multiple of N')
            end
            num_equilibria = N_alpha / N;
            if problem.nU > 0
                if size(u_star,1) ~= problem.nU || size(u_star,2) ~= num_equilibria
                    error('penalization of alpha_eq: wrong dimension on u_star')
                end
                u_star = renormalize_m1p1(u_star,normaliz.u_min,normaliz.u_max);
            end
            AlphaGrid = [];
                        
%             alpha_norm = ones(N_alpha,1);
%             if ~useG
%                 for i_point = 1:num_equilibria
%                     alpha_norm((i_point-1)*N+1:(i_point-1)*N+problem.nY) = (normaliz.y_max-normaliz.y_min)/2;
%                 end
%             end
            pen_handler.alpha_norm.max =  ones(N_alpha,1);
            pen_handler.alpha_norm.min = -ones(N_alpha,1);
            if ~useG
                for i_point = 1:num_equilibria
                    pen_handler.alpha_norm.max((i_point-1)*N+1:(i_point-1)*N+problem.nY) = normaliz.y_max;
                    pen_handler.alpha_norm.min((i_point-1)*N+1:(i_point-1)*N+problem.nY) = normaliz.y_min;
                end
            end
    
            pen_handler.size_pen = nS*N*num_equilibria;
            pen_handler.norm_pen = 1/pen_handler.size_pen;
            
%             pen_handler.get_alpha_to_alpha_handler = @get_alpha_to_alpha_handler;
        end
    end

    pen_handler.init_epoch = @init_epoch;
    function init_epoch(pen_handler,paramsF,paramsG,Alpha,IC)
        if pen_handler.pen_active
            AlphaGrid = Alpha;
        end
    end

    pen_handler.evaluate = @evaluate;
    function out = evaluate(paramsF,paramsG,Alpha,IC,compRes,compGrad,compJac,pen_handler,numIter,modelclass)
        out.dummy = 0;
        if pen_handler.coef > 0
            
            res = zeros(nS,num_equilibria,N);
            if compJac || compGrad
                jac = zeros(nS,num_equilibria,N,nwF);
            end
			for iS = 1:nS
                for i_point = 1:num_equilibria
                    x_curr = AlphaGrid((i_point-1)*N+1:i_point*N, iS);                    
                    res(iS,i_point,:) = modelclass.eval_f(x_curr,u_star(:, i_point),AlphaGrid(:, iS),paramsF,[]);
                    if compJac || compGrad
                        [jac(iS,i_point,:,:),~,~] = modelclass.eval_sensitivity_f();
                    end
                end
			end

            res = reshape(res,[],1);
            if compJac || compGrad
                jac = reshape(jac,[],nwF);
            end
            
            out.errors = sum(res.^2);
            if compGrad
                out.gradF{1} = 2*jac'*res;
            end
            if compRes
                out.res{1} = res;
            end
            if compJac
                out.jacF{1} = jac;
            end
        end          
    end

%     function hdl = get_alpha_to_alpha_handler(HFmod)
%         hdl = @(a) alpha_to_alpha(HFmod,a);
%     end
% 
%     function a_new = alpha_to_alpha(HFmod, a)
%         if N > 1 || useG
%             error('not (yet) implemented!')
%         end
%         x_grid = linspace(HFmod.problem.y_min, HFmod.problem.y_max, 1e2);
%         for i_point = 1:num_equilibria
%             a_new{i_point} = [];
%             f_grid = zeros(1,length(x_grid));
%             for iX = 1:length(x_grid)
%                 f_grid(iX) = HFmod.f_alpha(x_grid,u_star(:, i_point),a);
%                 if f_grid(iX) == 0
%                     a_new{i_point} = [a_new{i_point} x_grid(iX)];
%                 elseif iX > 1 && f_grid(iX)*f_grid(iX-1) < 0
%                     x_0 = (f_grid(iX)*x_grid(iX-1) - f_grid(iX-1)*x_grid(iX)) / (f_grid(iX) - f_grid(iX-1));
%                     a_new{i_point} = [a_new{i_point} x_0];
%                 end
%             end
%         end
%     end
end