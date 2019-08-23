function pen_handler = penalization_dfdalpha(optionsfile)

    dfdaPen_s           = iniread(optionsfile,'Penalizations','pen_dfda','s','0');
    dfdaPen_Npt         = iniread(optionsfile,'Penalizations','pen_dfda_Npt','d','10');
    
    pen_handler.num_pen = 1;
    pen_handler.names = {'ERRdfda'};
    pen_handler.coefficient_string = {dfdaPen_s};
    pen_handler.dependence_f = 1;
    pen_handler.dependence_g = 0;
    pen_handler.dependence_a = 0;
    pen_handler.dependence_i = 0;
    
    dfda_pt = [];
    AlphaGrid = [];
    N = [];
    N_alpha = [];
    nwF = [];
    nS = [];
       
    pen_handler.initialize = @initialize;
    function pen_handler = initialize(pen_handler,problem,N_,N_alpha_,useG,nS_,nwF_,nwG,nAL,nIC,misc)
        if pen_handler.pen_active
            N = N_;
            N_alpha = N_alpha_;
            nwF = nwF_;
            nS = nS_;
            if problem.nU > 0
                error('penalization of dfda works only if nU = 0')
            end
            rng('default');
            dfda_pt = 2*lhsdesign(nS*dfdaPen_Npt,N+problem.nU)-1; % *dfdaPen_Npt points in [-1,1]^(N+nU)
            dfda_pt = reshape(dfda_pt,nS,dfdaPen_Npt,[]);
            AlphaGrid = [];
    
            pen_handler.size_pen = nS*dfdaPen_Npt*N*N_alpha;
            pen_handler.norm_pen = 1/pen_handler.size_pen;
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
            ResDfda = zeros(nS,dfdaPen_Npt,N,N_alpha);
            if compJac || compGrad
                JacDfda = zeros(nS,dfdaPen_Npt,N,N_alpha,nwF);
            end
            for iS = 1:nS
                for iPt = 1:dfdaPen_Npt
                    [Da,Daw] = modelclass.eval_sensitivity_dfdaPen(paramsF,reshape(dfda_pt(iS,iPt,:),[],1),AlphaGrid(:,iS));
                    ResDfda(iS,iPt,:,:) = Da;
                    if compJac || compGrad
                        JacDfda(iS,iPt,:,:,:) = Daw;
                    end
                end
            end
            ResDfda = reshape(ResDfda,[],1);
            if compJac || compGrad
                JacDfda = reshape(JacDfda,[],nwF);
            end
            
            out.errors = sum(ResDfda.^2);
            if compGrad
                out.gradF{1} = 2*JacDfda'*ResDfda;
            end
            if compRes
                out.res{1} = ResDfda;
            end
            if compJac
                out.jacF{1} = JacDfda;
            end
        end          
    end

end