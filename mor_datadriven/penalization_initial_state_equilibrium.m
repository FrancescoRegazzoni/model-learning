function pen_handler = penalization_initial_state_equilibrium(optionsfile)

    ref_s = iniread(optionsfile,'Penalizations','pen_initial_state_eq','s','0');
    
    pen_handler.num_pen = 1;
    pen_handler.names = {'ERRinitStateEq'};
    pen_handler.coefficient_string = {ref_s};
    pen_handler.dependence_f = 1;
    pen_handler.dependence_g = 0;
    pen_handler.dependence_a = 1;
    pen_handler.dependence_i = 1;
       
    N = [];
    N_alpha = [];
    nwF = [];
    nAL = [];
    nIC = [];
    nS = [];
    u0_tab = [];
    
    pen_handler.initialize = @initialize;
    function pen_handler = initialize(pen_handler,problem,N_,N_alpha_,useG,nS_,nwF_,nwG,nAL_,nIC_,misc)
        if pen_handler.pen_active(1) 
            N = N_;
            N_alpha = N_alpha_;
            nwF = nwF_;
            nAL = nAL_;
            nIC = nIC_;
            nS = nS_;   
            pen_handler.size_pen(1) = N*nS;
            pen_handler.norm_pen(1) = 1/pen_handler.size_pen(1);
            
            if problem.nU > 0
                for iS = 1:nS
                    u0_tab(:,iS) = misc.tr{iS}.uu(:,1);
                end
            end
            
        end
    end

    pen_handler.evaluate = @evaluate;
    function out = evaluate(paramsF,paramsG,Alpha,IC,compRes,compGrad,compJac,pen_handler,numIter,modelclass)
        out.dummy = 0;
        if pen_handler.coef > 0
                        
            res = zeros(nS,N);
            if compJac || compGrad
                jacF = zeros(nS,N,nwF);
                jacA = zeros(nS,N,N_alpha,nS);
                jacI = zeros(nS,N,N,nS);
            end
            
            for iS = 1:nS
                res(iS,:) = modelclass.eval_f(IC(:,iS),u0_tab(:,iS),Alpha(:,iS),paramsF,1);
                if compJac || compGrad
                    [jacF(iS,:,:),jacI(iS,:,:,iS),jacA(iS,:,:,iS)] = modelclass.eval_sensitivity_f();
                end
            end

            res = reshape(res,[],1);
            if compJac || compGrad
                jacF = reshape(jacF,[],nwF);
                jacA = reshape(jacA,[],nAL);
                jacI = reshape(jacI,[],nIC);
            end
            
            out.errors = sum(res.^2);
            if compRes
                out.res{1} = res;
            end
            if compGrad
                out.gradF{1} = 2*jacF'*res;
                out.gradA{1} = 2*jacA'*res;
                out.gradI{1} = 2*jacI'*res;
            end
            if compJac
                out.jacF{1} = jacF;
                out.jacA{1} = jacA;
                out.jacI{1} = jacI;
            end   
        end
    end

end