function pen_handler = penalization_rhs_increasing(optionsfile)

    Pen_s           = iniread(optionsfile,'Penalizations','pen_rhs_incr','s','0');
    
    pen_handler.num_pen = 1;
    pen_handler.names = {'ERRrhsIncr'};
    pen_handler.coefficient_string = {Pen_s};
    pen_handler.dependence_f = 1;
    pen_handler.dependence_g = 0;
    pen_handler.dependence_a = 0;
    pen_handler.dependence_i = 0;
    
    N = [];
    nwF = [];
    N_alpha = [];
    AlphaBounds = [];
    nsamples = 1e2;
    
    X = [];
       
    pen_handler.initialize = @initialize;
    function pen_handler = initialize(pen_handler,problem,N_,N_alpha_,useG_,nS,nwF_,nwG,nAL,nIC,misc)
        if pen_handler.pen_active
            N = N_;
            nwF = nwF_;
            N_alpha = N_alpha_;
            if problem.nU ~=1 || N ~=1
                error('penalization_steady_state_increasing: wrong dimensions')
            end
            if useG_
                error('penalization_steady_state_increasing works only for useG_ = false')
            end
            
            pen_handler.size_pen = nsamples;
            pen_handler.norm_pen = 1/pen_handler.size_pen;
            
            X = lhsdesign(nsamples, N_alpha + 2);
        end
    end

    pen_handler.init_epoch = @init_epoch;
    function init_epoch(pen_handler,paramsF,paramsG,Alpha,IC)
        if pen_handler.pen_active
            AlphaBounds = [min(Alpha,[],2) max(Alpha,[],2)];
        end
    end

    pen_handler.evaluate = @evaluate;
    function out = evaluate(paramsF,paramsG,Alpha,IC,compRes,compGrad,compJac,pen_handler,numIter,modelclass)
        out.dummy = 0;
        if pen_handler.coef > 0
            
            ep = 1e-2;
%             penal_func = @(x) (x<=ep).*((x/ep-1).^2);
%             d_penal_func = @(x) 2/ep*(x<=ep).*(x/ep-1);
            penal_func = @(x) max(0,ep-x);
            d_penal_func = @(x) -(x<=ep);
            
            res = zeros(nsamples,1);
            if compJac || compGrad
                jac = zeros(nsamples,nwF);
            end
            
            for i_point = 1:nsamples
                x_curr = X(i_point, 1)*2 -1;
                u_curr = X(i_point, 2)*2 -1;
                a_curr = X(i_point, 3:end)'.*(AlphaBounds(:,2) - AlphaBounds(:,1)) + AlphaBounds(:,1);
                
                [~,Dx,~,Dxw] = modelclass.eval_sensitivity_steady_state_incr(paramsF,x_curr,u_curr,a_curr);
                
                res(i_point) = penal_func(Dx);
                if compJac || compGrad
                    jac(i_point,:) = d_penal_func(Dx) * Dxw;
                end
            end

%             res = reshape(res,[],1);
%             if compJac || compGrad
%                 jac = reshape(jac,[],nwF);
%             end
            
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
end