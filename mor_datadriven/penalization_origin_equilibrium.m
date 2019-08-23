function pen_handler = penalization_origin_equilibrium(optionsfile)

    origEqPen_s         = iniread(optionsfile,'Penalizations','pen_origstable_eq','s','0');
    origStabPen_s       = iniread(optionsfile,'Penalizations','pen_origstable_stab','s','0');
    origStabTh_s        = iniread(optionsfile,'Penalizations','pen_origstable_stab_th','s');
    origStabU           = iniread(optionsfile,'Penalizations','pen_origstable_u');
    origStabY           = iniread(optionsfile,'Penalizations','pen_origstable_y');
    
    pen_handler.num_pen = 2;
    pen_handler.names = {'ERRorigEq','ERRorigStab'};
    pen_handler.coefficient_string = {origEqPen_s, origStabPen_s};
    pen_handler.dependence_f = [1,1];
    pen_handler.dependence_g = [0,0];
    pen_handler.dependence_a = [0,0];
    pen_handler.dependence_i = [0,0];
    
    origStabTh_fun  = inline(origStabTh_s ,'k');
   
    nP = [];
    nResOrig = [];
    N = [];
    nwF = [];
    nY = [];
    xhat_end_curr = [];
    
    pen_handler.initialize = @initialize;
    function pen_handler = initialize(pen_handler,problem,N_,N_alpha,useG,nS,nwF_,nwG,nAL,nIC,misc)
        if pen_handler.pen_active(1) || pen_handler.pen_active(2)
            N = N_;
            nwF = nwF_;
            nY = problem.nY;
            if N_alpha > 0
                error('OriginStable penalization works only if N_alpha = 0')
            end                
            nP = size(origStabU,2);
            nResOrig = nP*N;
            pen_handler.norm_pen(1) = 1/nResOrig;
            pen_handler.norm_pen(2) = 1/nResOrig;
            pen_handler.size_pen(1) = nResOrig;
            pen_handler.size_pen(2) = nResOrig;
            
            if problem.nU > 0
                origStabU = renormalize_m1p1(origStabU,misc.normaliz.u_min,misc.normaliz.u_max);
            end
            origStabY = renormalize_m1p1(origStabY,misc.normaliz.y_min,misc.normaliz.y_max);

            for iP = 1:nP
                xhat_end_curr(:,iP) = misc.xhat_end; 
                if ~useG
                    xhat_end_curr(1:nY,iP) = origStabY(:,iP);
                end
            end
        end
    end

    pen_handler.evaluate = @evaluate;
    function out = evaluate(paramsF,paramsG,Alpha,IC,compRes,compGrad,compJac,pen_handler,numIter,modelclass)
        out.dummy = 0;        
        if pen_handler.coef(1) > 0 || pen_handler.coef(2) > 0
            origStabTh  = origStabTh_fun(numIter); 
    %         PhiStab  = @(w) (w>origStabTh).*(w/origStabTh - 1).^2;
    %         dPhiStab = @(w) (w>origStabTh).*(w/origStabTh - 1).*(2/origStabTh);   
            PhiStab  = @(w) (w>origStabTh).*(1 - w/origStabTh);
            dPhiStab = @(w) (w>origStabTh)/(-origStabTh);        

            ResOrigEq = zeros(nResOrig,1);
            ResOrigStab = zeros(nResOrig,1);
            if compJac || compGrad
                JacOrigEqF = zeros(nResOrig,nwF);
                JacOrigStabF = zeros(nResOrig,nwF);
            end

    %         xOrig = xhat_init;
    %         if ~useG
    %             xOrig(1:nY) = Yref;
    %         end

            for iP = 1:nP
                [outputF,Dx,Dw,Dwx] = modelclass.eval_sensitivity_origStabPen(paramsF,xhat_end_curr(:,iP),origStabU(:,iP));
                idxCurr = N*(iP-1)+1:N*iP;

                if pen_handler.coef(1) > 0
                    ResOrigEq(idxCurr) = outputF;
                end
                if pen_handler.coef(2) > 0
                    Jorig = Dx;
                    [V,D,W] = eig(Jorig);
                    lam = diag(D);
    %                 ResOrigStab(idxCurr) = PhiStab(diag(Jorig));
                    ResOrigStab(idxCurr) = PhiStab(real(lam));
                end

                if compJac || compGrad
                    if pen_handler.coef(1) > 0
                        JacOrigEqF(idxCurr,:) = Dw;
                    end
                    if pen_handler.coef(2) > 0
    %                     for iX = 1:N
    %                         JacOrigStabF(N*(iP-1)+iX,:) = dPhiStab(doutdinput(iX,nU+iX))* ...
    %                             [permute(ddoutdinputdw(iX,nU+iX,:),[1 3 2]) permute(ddoutdinputdtheta(iX,nU+iX,:),[1 3 2])];
    %                     end
                        for iDesign = 1:nwF
                            E = Dwx(:,nU+1:end,iDesign);
                            for iX = 1:N
                                JacOrigStabF(N*(iP-1)+iX,iDesign) = dPhiStab(real(lam(iX)))* ...
                                    real((W(:,iX)'*E*V(:,iX))/(W(:,iX)'*V(:,iX)));
                            end
                        end
                    end
                end
            end

            out.errors(1) = sum(ResOrigEq.^2);
            out.errors(2) = sum(ResOrigStab.^2);
            if compGrad
                out.gradF{1} = 2*JacOrigEqF'*ResOrigEq;
                out.gradF{2} = 2*JacOrigStabF'*ResOrigStab;
            end
            if compRes
                out.res{1} = ResOrigEq;
                out.res{2} = ResOrigStab;
            end
            if compJac
                out.jacF{1} = JacOrigEqF;
                out.jacF{2} = JacOrigStabF;
            end     
        end          
    end

end