function pen_handler = penalization_ref_equilibrium(optionsfile)

    refEqPen_s         = iniread(optionsfile,'Penalizations','pen_ref_eq','s','0');
    
    pen_handler.num_pen = 1;
    pen_handler.names = {'ERRrefEq'};
    pen_handler.coefficient_string = {refEqPen_s};
    pen_handler.dependence_f = 1;
    pen_handler.dependence_g = 0;
    pen_handler.dependence_a = 0;
    pen_handler.dependence_i = 0;
       
    N = [];
    nwF = [];
    u_ref = [];
    x_ref = [];
    
    pen_handler.initialize = @initialize;
    function pen_handler = initialize(pen_handler,problem,N_,N_alpha,useG,nS,nwF_,nwG,nAL,nIC,misc)
        if pen_handler.pen_active(1) 
            N = N_;
            nwF = nwF_;
            if N_alpha > 0
                error('RefEq penalization works only if N_alpha = 0')
            end                
            pen_handler.norm_pen(1) = 1/N;
            pen_handler.size_pen(1) = N;
            
            if problem.nU > 0
                u_ref = renormalize_m1p1(problem.u_ref,misc.normaliz.u_min,misc.normaliz.u_max);
            end
            x_ref = misc.xhat_end;
        end
    end

    pen_handler.evaluate = @evaluate;
    function out = evaluate(paramsF,paramsG,Alpha,IC,compRes,compGrad,compJac,pen_handler,numIter,modelclass)
        out.dummy = 0;        
        if pen_handler.coef(1) > 0
            
            ResRefEq = zeros(N,1);
            if compJac || compGrad
                JacRefEqF = zeros(N,nwF);
            end

            [outputF,~,Dw,~] = modelclass.eval_sensitivity_origStabPen(paramsF,x_ref,u_ref);

            if pen_handler.coef(1) > 0
                ResRefEq = outputF;
            end

            if compJac || compGrad
                if pen_handler.coef(1) > 0
                    JacRefEqF = Dw;
                end
            end

            out.errors(1) = sum(ResRefEq.^2);
            if compGrad
                out.gradF{1} = 2*JacRefEqF'*ResRefEq;
            end
            if compRes
                out.res{1} = ResRefEq;
            end
            if compJac
                out.jacF{1} = JacRefEqF;
            end     
        end          
    end

end