	function pen_handler = penalization_f_constrained_wrt_alpha(optionsfile)

    penal_s           = iniread(optionsfile,'Penalizations','pen_fconstr','s','0');
    x_star            = iniread(optionsfile,'Penalizations','pen_fconstr_xstar','d','0');
    u_star            = iniread(optionsfile,'Penalizations','pen_fconstr_ustar','d','0');
    Npt_alpha         = iniread(optionsfile,'Penalizations','pen_fconstr_Npt','d','10');
    grid_border_width = iniread(optionsfile,'Penalizations','pen_fconstr_grid_bd_wd','d','.2');
    initalphadir      = iniread(optionsfile,'Penalizations','pen_fconstr_initalphadir','s','');
    
    pen_handler.num_pen = 1;
    pen_handler.names = {'ERRfconstr'};
    pen_handler.coefficient_string = {penal_s};
    pen_handler.dependence_f = 1;
    pen_handler.dependence_g = 0;
    pen_handler.dependence_a = 0;
    pen_handler.dependence_i = 0;
	
	AlphaGrid = [];
	N_alpha = [];
    nwF = [];
    x_star_orig = [];
    u_star_orig = [];
%     normaliz = [];
    alpha_norm = [];
       
    pen_handler.initialize = @initialize;
    function pen_handler = initialize(pen_handler,problem,N,N_alpha_,useG,nS,nwF_,nwG,nAL,nIC,misc)
        if pen_handler.pen_active
			N_alpha = N_alpha_;
            nwF = nwF_;
            normaliz = misc.normaliz;
            if N ~= N_alpha
                error('This penalization works only if N == N_alpha')
            end
            u_star = adapt_dimension(u_star,problem.nU);
            x_star = adapt_dimension(x_star,N);
            
            x_star_orig = x_star;
            u_star_orig = u_star;
            
            if problem.nU > 0
                u_star = renormalize_m1p1(u_star,normaliz.u_min,normaliz.u_max);
            end
            if ~useG
                x_star(1:problem.nY) = renormalize_m1p1(x_star(1:problem.nY),normaliz.y_min,normaliz.y_max);
            end
            
            pen_handler.size_pen = Npt_alpha*N_alpha;
            pen_handler.norm_pen = 1/pen_handler.size_pen;
            
            alpha_norm = ones(N,1);
            if ~useG
                alpha_norm(1:problem.nY) = (normaliz.y_max-normaliz.y_min)/2;
            end
            alpha_norm = alpha_norm./normaliz.t_norm;
            pen_handler.alpha_norm = alpha_norm;
            
            pen_handler.get_alpha_to_alpha_handler = @get_alpha_to_alpha_handler;
        end
    end

    pen_handler.fix_alpha = @fix_alpha;
    function Alpha_0 = fix_alpha(pen_handler,problem,Alpha_0)
        if pen_handler.pen_active && ~isempty(initalphadir)
            ANNmod = read_model_fromfile(problem,initalphadir);
            for iS = 1:size(Alpha_0,2)
                Alpha_0(:,iS) = ANNmod.f_alpha(x_star_orig,u_star_orig,ANNmod.alpha_learned(iS))./alpha_norm;
            end
        end
    end
	
    pen_handler.init_epoch = @init_epoch;
    function init_epoch(pen_handler,paramsF,paramsG,Alpha,IC)
        if pen_handler.pen_active
            a_max = max(Alpha,[],2)';
            a_min = min(Alpha,[],2)';
            a_max_ext = a_max + max(.1,grid_border_width*(a_max-a_min));
            a_min_ext = a_min - max(.1,grid_border_width*(a_max-a_min));
            rng('default')
            AlphaGrid = a_min_ext + (a_max_ext-a_min_ext).*lhsdesign(Npt_alpha,N_alpha);
        end
    end

    pen_handler.evaluate = @evaluate;
    function out = evaluate(paramsF,paramsG,Alpha,IC,compRes,compGrad,compJac,pen_handler,numIter,modelclass)
        out.dummy = 0;
        if pen_handler.coef > 0
            res = zeros(Npt_alpha,N_alpha);
            if compJac || compGrad
                jac = zeros(Npt_alpha,N_alpha,nwF);
            end
			for iP = 1:Npt_alpha
				a = AlphaGrid(iP,:)';
				res(iP,:) = (modelclass.eval_f(x_star,u_star,a,paramsF,1) - a)'; 
				if compJac || compGrad
					[jac(iP,:,:),~,~] = modelclass.eval_sensitivity_f();
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
    
    function hdl = get_alpha_to_alpha_handler(HFmod)
        hdl = @(a) HFmod.f_alpha(x_star_orig,u_star_orig,a);
    end

end