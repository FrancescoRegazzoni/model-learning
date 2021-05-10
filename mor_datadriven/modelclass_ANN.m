function modelclass = modelclass_ANN(optionsfile,problem,N,N_alpha,useG,constrain_type)

    % constrain_type:
    % 0: free; 
    % 1: f(x^*,u^*,a) = a
    % 2: f(x^*_1,u^*_1,a) = [I,0]a, f(x^*_2,u^*_2,a) = [0,I]a
    % 3: f(x_ref,u_ref) = 0
    % 4: f(a_1,u^*_1,a) = 0, f(a_2,u^*_2,a) = 0
    
    
    if nargin == 5
        constrain_type = 0; 
    end
    
    constrained_alpha = 0;
    constrained_alpha_single = 0;
    constrained_alpha_double = 0;
    constrained_ref = 0;
    constrained_steadystate = 0;
    switch constrain_type
        case 0 % free
            constrained = 0;
        case 1 % f(x^*,u^*,a) = a
            constrained = 1;
            constrained_alpha = 1;
            constrained_alpha_single = 1;
        case 2 % f(x^*_1,u^*_1,a) = [I,0]a, f(x^*_2,u^*_2,a) = [0,I]a
            constrained = 1;
            constrained_alpha = 1;
            constrained_alpha_double = 1;
        case 3 % f(x_ref,u_ref) = 0
            constrained = 1;
            constrained_ref = 1;
        case 4 % f(a_1,u^*_1,a) = 0, f(a_2,u^*_2,a) = 0
            constrained = 1;
            constrained_steadystate = 1;
    end
    
    layF                = iniread(optionsfile,'ANN','layF');
    layG                = iniread(optionsfile,'ANN','layG');
    BetaOutput          = iniread(optionsfile,'ANN','BetaOutput','d',1);
    init_type           = iniread(optionsfile,'ANN','init_type','s','shuffle');
    init_file           = iniread(optionsfile,'ANN','init_file','s');
    init_iter           = iniread(optionsfile,'ANN','init_iter','s','end');
    
    dynamics_type       = iniread(optionsfile,'Model','dynamics_type','s','internal');

    switch dynamics_type
        case 'internal'
            int_dyn = 1;
        case 'external'
            int_dyn = 0;
        otherwise
            error('dynamics type not recognised')
    end
    
    baseopt = get_base_options();
    BaseDir = baseopt.BaseDir;
    SaveDir = strcat(BaseDir,'/',problem.dir_nets,'/');
    
    nU = problem.nU;
    nY = problem.nY;    
    
    if constrained_alpha
        initalphadir      = iniread(optionsfile,'Penalizations','pen_fconstr_initalphadir','s','');
    
        x_star            = iniread(optionsfile,'Penalizations','pen_fconstr_xstar','d','0');
        u_star            = iniread(optionsfile,'Penalizations','pen_fconstr_ustar','d','0');
        if problem.nU > 1 && size(u_star,1) == 1
            u_star = u_star*ones(problem.nU,1);
        elseif problem.nU == 0
            u_star = [];
        end
        if N > 1 && size(x_star,1) == 1
            x_star = x_star*ones(N,1);
        end
    end
    if constrained_alpha_double
        x_star_2          = iniread(optionsfile,'Penalizations','pen_fconstr_xstar_2','d','0');
        u_star_2          = iniread(optionsfile,'Penalizations','pen_fconstr_ustar_2','d','0');

        if problem.nU > 1 && size(u_star_2,1) == 1
            u_star_2 = u_star_2*ones(problem.nU,1);
        elseif problem.nU == 0
            u_star_2 = [];
        end
        if N > 1 && size(x_star_2,1) == 1
            x_star_2 = x_star_2*ones(N,1);
        end
    end
    if constrained_ref
        x_star = [];
        modelclass.set_xref = @set_xref;
        u_star = problem.u_ref;
    end
    if constrained_steadystate
        u_star          = iniread(optionsfile,'Penalizations','pen_alpha_eq_ustar','d','0');
    end
    function set_xref(x_reference)
        x_star = x_reference;
    end
    
    %% Check compatibility
    if ~int_dyn
%         if alphaPen_act || alphaRaise_act || ratioPen_act 
%             error('penalizations incompatible with external dynamics')
%         end
        if mod(N,nY) ~= 0
            error('with the external dynamics approach, N must be multiple of nY')
        end
        if useG
            error('useG option incompatible with external dynamics')
        end
        mod_ord = N/nY; % model order
    end
    
    if constrained_alpha_single && N ~= N_alpha
        error('N is different from N_alpha')
    end
    if constrained_alpha_double && 2*N ~= N_alpha
        error('2N is different from N_alpha')
    end
    if constrained_alpha && ~int_dyn
        error('alpha constrain requires int_dyn')
    end
    if constrained_ref && N_alpha > 0
        error('constrained_ref is available only if N_alpha = 0')
    end
    if constrained_alpha_double
        if norm(x_star - x_star_2) + norm(u_star - u_star_2) < 1e-8
            error('the two reference point are too close')
        end
    end
    if constrained_steadystate
        if N ~= 1 || N_alpha ~= 2 || problem.nU ~= 1
            error('wrong dimensions for model class.')
        end
        if useG
            error('constrained_steadystate incompatible with useG.')
        end
        if abs(u_star(1) - u_star(2)) < 1e-8
            error('the two reference inputs are too close')
        end
    end

    %% Neural Network initialization

    f = @(x) tanh(x);
    df = @(x) sech(x).^2;
    ddf = @(x) 2*tanh(x).*(tanh(x).^2 - 1);

    % f = @(x) 2./(1+exp(-x))-1;
    % df = @(x) 2*exp(-x)./(1+exp(-x)).^2;

    % Tf=.5;
    % f2 = @(x) 2./(1+exp(-x/Tf))-1;
    % df2 = @(x) 2/Tf*exp(-x/Tf)./(1+exp(-x/Tf)).^2;

    %NewInit = strcmp(init_file,'') || size(init_file,2)==0;
    if int_dyn
        numnF = [nU+N+N_alpha;layF;N];
    else
        numnF = [nU+N+N_alpha;layF;nY];
    end
    if useG
        numnG = [N;layG;nY];
    end
    if ~strcmp(init_type,'file') && ~strcmp(init_type,'f')
        switch(init_type)
            case {'shuffle','s'}
                rng('shuffle')
            case {'default','d'}
                rng('default')
            otherwise
                rng(str2double(init_type))
        end
        [wF0,thetaF0] = BuildRandomANN(numnF);
        if useG
            [wG0,thetaG0] = BuildRandomANN(numnG);
        end

    %     lastlayT = length(thetaF0)-numnF(end)+1:length(thetaF0);
    %     lastlayW = length(wF0)-numnF(end)*numnF(end-1)+1:length(wF0);
    %     thetaF0(lastlayT) = 0;
    %     wF0(lastlayW) = wF0(lastlayW)*1e-2;

    else
        init_data_file = strcat(SaveDir,init_file,'/network.mat');
        oldstylesave = 1;
        if ~exist(init_data_file,'file')
            oldstylesave = 0;
            init_data_file = strcat(SaveDir,init_file,'/data.mat');
        end        
        init_data = load(init_data_file,'numnF','BetaOutput');
        if ~isequal(numnF,init_data.numnF) || ~isequal(BetaOutput,init_data.BetaOutput)
            error('Inconsistent initialization')
        end
        if useG init_dataG = load(init_data_file,'numnG'); end
        if oldstylesave
            if strcmp(init_iter,'end')
                load(init_data_file,'wF','thetaF');
                wF0 = wF; thetaF0 = thetaF; 
                if useG
                    load(init_data_file,'wG','thetaG');
                    wG0 = wG; thetaG0 = thetaG; 
                end
            elseif strcmp(init_iter,'0')
                load(init_data_file,'wF0','thetaF0');
                if useG
                    load(init_data_file,'wG0','thetaG0');
                end
            else
                init_iterVal = str2num(init_iter);
                load(init_data_file,'wF_hist','thetaF_hist');
                wF0 = wF_hist(:,init_iterVal); 
                thetaF0 = thetaF_hist(:,init_iterVal); 
                if useG
                load(init_data_file,'wG_hist','thetaG_hist');
                    wG0 = wG_hist(:,init_iterVal); 
                    thetaG0 = thetaG_hist(:,init_iterVal); 
                end
            end
        else
            if strcmp(init_iter,'end')
                init_net_file = strcat(SaveDir,init_file,'/net.mat');
            elseif strcmp(init_iter,'0')
                init_net_file = strcat(SaveDir,init_file,'/networks/net_000000.mat');
            else
                init_iterVal = str2num(init_iter);
                init_net_file = strcat(SaveDir,init_file,'/networks/',sprintf('net_%06d.mat',init_iterVal));
            end    
            init_net = load(init_net_file,'wF','thetaF');
            wF0 = init_net.wF; thetaF0 = init_net.thetaF; 
            if useG
                init_net = load(init_net_file,'wG','thetaG');
                wG0 = init_net.wG; thetaG0 = init_net.thetaG; 
            end
        end
    end

    nwFw = length(wF0);
    nwFtheta = length(thetaF0);
    nwF = nwFw + nwFtheta;
    if useG
        nwGw = length(wG0);
        nwGtheta = length(thetaG0);
        nwG = nwGw + nwGtheta;
    end

    modelclass.nwF = nwF;
    modelclass.paramsF0 = [wF0;thetaF0];
    if constrained
        modelclass.nwF = modelclass.nwF - N;
        modelclass.paramsF0 = modelclass.paramsF0(1:end-N);
    end
    if useG
        modelclass.nwG = nwG;
        modelclass.paramsG0 = [wG0;thetaG0];
    end

    %% Labels
    if int_dyn
        modelclass.name_modelclass = '_int';
    else
        modelclass.name_modelclass = '_ext';
    end    
    modelclass.laystrF = sprintf('%d-',layF);
    modelclass.laystrF = ['_hlayF' modelclass.laystrF(1:end-1)];
    if useG
        modelclass.laystrG = sprintf('%d-',layG);
        modelclass.laystrG = ['_hlayG' modelclass.laystrG(1:end-1)];
    end
    
    %% inizial normalization
    if constrained_alpha
        modelclass.apply_normalization = @apply_normalization_alpha_constrain;
    end    
    if constrained_alpha
        x_star_orig = x_star;
        u_star_orig = u_star;
    end
    if constrained_alpha_double
        x_star_orig_2 = x_star_2;
        u_star_orig_2 = u_star_2;
    end
    alpha_norm = [];
    function alpha_norm_ret = apply_normalization_alpha_constrain(nr)
        if nU > 0
%             u_star = nr.u_norm(u_star);
            u_star = renormalize_m1p1(u_star,nr.u_min,nr.u_max);
        end
        if ~useG
%                 x_star(1:nY) = nr.y_norm(x_star(1:nY));
            x_star(1:nY) = renormalize_m1p1(x_star(1:nY),nr.y_min,nr.y_max);
            if constrained_alpha_double
%                 x_star_2(1:nY) = nr.y_norm(x_star_2(1:nY));
                x_star_2(1:nY) = renormalize_m1p1(x_star_2(1:nY),nr.y_min,nr.y_max);
            end
        end       
        alpha_norm_ret = ones(N_alpha,1);
        if ~useG
            alpha_norm_ret(1:nY) = (nr.y_max-nr.y_min)/2;
            if constrained_alpha_double
                alpha_norm_ret(N+1:2*nY) = (nr.y_max-nr.y_min)/2;
            end
        end
        alpha_norm_ret = alpha_norm_ret./nr.t_norm; 
        alpha_norm = alpha_norm_ret;
    end

    if constrained_ref
        modelclass.apply_normalization = @apply_normalization_ref_constrain;
    end
    function alpha_norm_ret = apply_normalization_ref_constrain(nr)
        if nU > 0
            u_star = renormalize_m1p1(u_star,nr.u_min,nr.u_max);
        end
        alpha_norm_ret = 0;
    end

    
    if constrained_steadystate
        modelclass.apply_normalization = @apply_normalization_steadystate;
    end
    function alpha_norm_ret = apply_normalization_steadystate(nr)
        u_star = renormalize_m1p1(u_star,nr.u_min,nr.u_max);
        alpha_norm_ret.min = ones(2,1) * nr.y_min;
        alpha_norm_ret.max = ones(2,1) * nr.y_max;
    end
        
    
    %% Renormalization
    if int_dyn
        if useG
            modelclass.idx_hidden_states = 1:N;
        else
            modelclass.idx_hidden_states = nY+1:N;
        end
    else
        modelclass.idx_hidden_states = [];
    end
    
    modelclass.renormalize_allows = 1;
    modelclass.renormalize = @renormalize;
    function [paramsF,paramsG] = renormalize(normFact,iX,paramsF,paramsG)
        if constrained
            paramsF = [paramsF;zeros(N,1)];
        end
        wF = paramsF(1:nwFw);
        thetaF = paramsF(nwFw+1:nwFw+nwFtheta);
        [wF,thetaF] = ANNmod_affinetransf_input(numnF,wF,thetaF,normFact,0,nU+iX);
        [wF,thetaF] = ANNmod_affinetransf_output(numnF,wF,thetaF,normFact,0,iX);
        if useG
            wG = paramsG(1:nwGw);
            thetaG = paramsG(nwGw+1:nwGw+nwGtheta);
            [wG,thetaG] = ANNmod_affinetransf_input(numnG,wG,thetaG,normFact,0,iX);
        end 
        if constrained
            paramsF = paramsF(1:end-N);
        end
    end

    if constrained_alpha
        modelclass.renormalize_allows = 0;
    end
        
    
    %% Plot options
    modelclass.plot_x = int_dyn;
    modelclass.x_idx_plot = modelclass.idx_hidden_states;
    
    %% Save functions
    modelclass.save_to_file_init = @save_to_file_init;
    function save_to_file_init(OutputFile)
        save(OutputFile,'f','df','numnF','BetaOutput',...
            'init_type','init_file','init_iter','int_dyn','',...
            '-append');
        if useG
            save(OutputFile,'numnG','-append');
        end
    end

    modelclass.save_to_file_paramsF = @save_to_file_paramsF;
    function save_to_file_paramsF(fileNet,paramsF)
        if constrained
            paramsF = [paramsF;zeros(N,1)];
        end
        wF = paramsF(1:nwFw);
        thetaF = paramsF(nwFw+1:nwFw+nwFtheta);
        save(fileNet,'thetaF','wF','-append');
    end
    
    modelclass.save_to_file_paramsG = @save_to_file_paramsG;
    function save_to_file_paramsG(fileNet,paramsG)
        wG = paramsG(1:nwGw);
        thetaG = paramsG(nwGw+1:nwGw+nwGtheta);
        save(fileNet,'thetaG','wG','-append');
    end

    %% Initial conditions (optional)   
    modelclass.fix_IC = @fix_IC;
    function IC_0 = fix_IC(IC_0,tr)
        if ~int_dyn
            for iS = 1:length(tr)
    %             IC_0(:,iS) = repmat(tr{iS}.yy(:,1),mod_ord,1);
                for iOrd = 1:mod_ord
                    IC_0((iOrd-1)*nY+1:iOrd*nY,iS) = tr{iS}.yy(:,1) - (iOrd-1)*(tr{iS}.yy(:,2)-tr{iS}.yy(:,1));
                end
            end 
        end
    end

    %% Model definition

    %modelclass.eval_f = @(x,u,a,w) EvaluateANN(numnF,w(1:nwFw),w(nwFw+1:nwFw+nwFtheta),[u;x;a],f,BetaOutput);
    

    % NB: we store precomputed values, because we are assuming that the
    % sentivity evaluation is called just after the corresponding evaluation
    curr_dt = [];
    curr_paramsF = [];
    curr_alphaF = [];
    curr_betaF = [];  
    curr_paramsG = []; 
    curr_alphaG = [];
    curr_betaG = [];   
    
    if constrained_alpha || constrained_ref
        curr_alphaF_base = [];
        curr_betaF_base = [];  
    end
    if constrained_alpha_double >= 2
        curr_alphaF_base_2 = [];
        curr_betaF_base_2 = [];  
        lambda = [];
        diff12 = [];
        Proj1 = [eye(N) zeros(N)];
        Proj2 = [zeros(N) eye(N)];
        vec_times_lambda = [];
    end
    if constrained_steadystate   
        outputF_base_1 = [];
        outputF_base_2 = [];
        curr_alphaF_base_1 = [];
        curr_betaF_base_1 = []; 
        curr_alphaF_base_2 = [];
        curr_betaF_base_2 = [];   
        lambda = [];
        d_lambda_dx = [];
        d_lambda_da1 = [];
        d_lambda_da2 = [];
    end
    modelclass.eval_f = @eval_f; 
    function outputF = eval_f(x,u,a,w,dt)
        if constrained
            w = [w;zeros(N,1)];
        end
        curr_dt = dt;
        curr_paramsF = w;
        [outputF,curr_alphaF,curr_betaF] = EvaluateANN(numnF,w(1:nwFw),w(nwFw+1:nwFw+nwFtheta),[u;x;a],f,BetaOutput);
        
        if constrained_alpha || constrained_ref
            [outputF_base  ,curr_alphaF_base  ,curr_betaF_base  ] = EvaluateANN(numnF,w(1:nwFw),w(nwFw+1:nwFw+nwFtheta),[u_star  ;x_star  ;a],f,BetaOutput);
        end
        if constrained_alpha_double
            [outputF_base_2,curr_alphaF_base_2,curr_betaF_base_2] = EvaluateANN(numnF,w(1:nwFw),w(nwFw+1:nwFw+nwFtheta),[u_star_2;x_star_2;a],f,BetaOutput);
        end
        if constrained_steadystate
            [outputF_base_1 ,curr_alphaF_base_1  ,curr_betaF_base_1 ] = EvaluateANN(numnF,w(1:nwFw),w(nwFw+1:nwFw+nwFtheta),[u_star(:,1)  ;a(1)  ;a],f,BetaOutput);
            [outputF_base_2 ,curr_alphaF_base_2  ,curr_betaF_base_2 ] = EvaluateANN(numnF,w(1:nwFw),w(nwFw+1:nwFw+nwFtheta),[u_star(:,2)  ;a(2)  ;a],f,BetaOutput);
        end
        
        if constrained_ref
            outputF = outputF - outputF_base;
        elseif constrained_alpha_single
            outputF = outputF - outputF_base + a;
        elseif constrained_alpha_double
            diff12 = [u_star_2;x_star_2] - [u_star;x_star];
            lambda = ([u;x]'*diff12) / (diff12'*diff12);
            vec_times_lambda = Proj2*a - Proj1*a + outputF_base - outputF_base_2;
            outputF = outputF + (Proj1*a - outputF_base) + lambda * vec_times_lambda;
        elseif constrained_steadystate
            denominator = ((a(2)-a(1))^2 + (u_star(:,2) - u_star(:,1))^2);
            lambda = ((x-a(1))*(a(2)-a(1)) + (u - u_star(:,1))*(u_star(:,2) - u_star(:,1)))/denominator;
            outputF = outputF - outputF_base_1 - (outputF_base_2-outputF_base_1) * lambda;
            d_lambda_dx = (a(2)-a(1))/denominator;
            d_lambda_da1 = (denominator*(2*a(1)-x-a(2)) +2*(a(2)-a(1))*((x-a(1))*(a(2)-a(1)) + (u - u_star(:,1))*(u_star(:,2) - u_star(:,1))))/denominator^2;
            d_lambda_da2 = (denominator*( x-a(1))       -2*(a(2)-a(1))*((x-a(1))*(a(2)-a(1)) + (u - u_star(:,1))*(u_star(:,2) - u_star(:,1))))/denominator^2;
        end
        
        if ~int_dyn
            outputF = [outputF; (x(1:end-nY,1)-x(nY+1:end,1))/dt];
        end
         
    end

    modelclass.eval_g = @eval_g; 
    function outputG = eval_g(x,w)
        curr_paramsG = w;
        [outputG,curr_alphaG,curr_betaG] = EvaluateANN(numnG,w(1:nwGw),w(nwGw+1:nwGw+nwGtheta),x,f,BetaOutput);
    end

    modelclass.eval_sensitivity_f = @eval_sensitivity_f;
    function [Dw,Dx,Da] = eval_sensitivity_f()            
        [doutdalphaF,doutdwF,doutdthetaF] = ANNSensitivityBP(numnF,curr_paramsF(1:nwFw),df,curr_alphaF,curr_betaF,BetaOutput);
        
        if constrained
            doutdthetaF = doutdthetaF(:,1:end-N);
        end
        if constrained_alpha || constrained_ref            
            [doutdalphaF_base  ,doutdwF_base  ,doutdthetaF_base  ] = ANNSensitivityBP(numnF,curr_paramsF(1:nwFw),df,curr_alphaF_base  ,curr_betaF_base  ,BetaOutput);
            doutdthetaF_base = doutdthetaF_base(:,1:end-N);
        end
        if constrained_alpha_double
            [doutdalphaF_base_2,doutdwF_base_2,doutdthetaF_base_2] = ANNSensitivityBP(numnF,curr_paramsF(1:nwFw),df,curr_alphaF_base_2,curr_betaF_base_2,BetaOutput);
            doutdthetaF_base_2 = doutdthetaF_base_2(:,1:end-N);
        end
        if constrained_steadystate
            [doutdalphaF_base_1,doutdwF_base_1,doutdthetaF_base_1] = ANNSensitivityBP(numnF,curr_paramsF(1:nwFw),df,curr_alphaF_base_1,curr_betaF_base_1,BetaOutput);
            doutdthetaF_base_1 = doutdthetaF_base_1(:,1:end-N);
            [doutdalphaF_base_2,doutdwF_base_2,doutdthetaF_base_2] = ANNSensitivityBP(numnF,curr_paramsF(1:nwFw),df,curr_alphaF_base_2,curr_betaF_base_2,BetaOutput);
            doutdthetaF_base_2 = doutdthetaF_base_2(:,1:end-N);
        end
        
        Da = [];
        if int_dyn
            Dw = [doutdwF doutdthetaF];
            Dx = doutdalphaF(:,nU+1:nU+N);
            if N_alpha > 0
                Da = doutdalphaF(:,nU+N+1:nU+N+N_alpha);
            end
            if constrained_ref
                Dw = Dw - [doutdwF_base doutdthetaF_base];
            end
            if constrained_alpha_single
                Dw = Dw - [doutdwF_base doutdthetaF_base];
                if N_alpha > 0
                    Da = Da - doutdalphaF_base(:,nU+N+1:nU+N+N_alpha) + eye(N);
                end
            end
            if constrained_alpha_double                
                Dw = Dw - (1-lambda)*[doutdwF_base doutdthetaF_base] - lambda*[doutdwF_base_2 doutdthetaF_base_2];
                Dx = Dx + vec_times_lambda * diff12(nU+1:end)' / (diff12'*diff12);
                if N_alpha > 0
                    Da = Da - (1-lambda)*doutdalphaF_base(:,nU+N+1:nU+N+N_alpha) - lambda*doutdalphaF_base_2(:,nU+N+1:nU+N+N_alpha) ...
                        + (1-lambda)*Proj1 + lambda*Proj2;
                end
            end
            if constrained_steadystate
                Dw = Dw - (1-lambda)*[doutdwF_base_1 doutdthetaF_base_1] - lambda*[doutdwF_base_2 doutdthetaF_base_2];
                Dx = Dx - (outputF_base_2-outputF_base_1) * d_lambda_dx;
                if N_alpha > 0
                    Da(:,1) = Da(:,1) - (1-lambda)*(doutdalphaF_base_1(:,nU+1:nU+N)+doutdalphaF_base_1(:,nU+N+1)) ...
                            -    lambda *doutdalphaF_base_2(:,nU+N+1) - (outputF_base_2-outputF_base_1)*d_lambda_da1;
                    Da(:,2) = Da(:,2) -    lambda *(doutdalphaF_base_2(:,nU+1:nU+N)+doutdalphaF_base_1(:,nU+N+2)) ...
                            - (1-lambda)*doutdalphaF_base_1(:,nU+N+2) - (outputF_base_2-outputF_base_1)*d_lambda_da2;
                end
            end
        else
            Dw(1:nY,:) = [doutdwF doutdthetaF];
            Dw(nY+1:N,:) = 0;
            Dx(1:nY,:) = doutdalphaF(:,nU+1:nU+N);
            for iOrd = 1:mod_ord-1
                Dx(iOrd*nY+1:(iOrd+1)*nY,(iOrd-1)*nY+1:(iOrd+1)*nY) = [eye(nY) -eye(nY)]/curr_dt;
            end
            if N_alpha > 0
                Da(1:nY,:) = doutdalphaF(:,nU+N+1:nU+N+N_alpha);
                Da(nY+1:N,:) = 0;
            end
        end
    end

    modelclass.eval_sensitivity_g = @eval_sensitivity_g;
    function [Dw,Dx] = eval_sensitivity_g()   
        [doutdalphaG,doutdwG,doutdthetaG] = ANNSensitivityBP(numnG,curr_paramsG(1:nwGw),df,curr_alphaG,curr_betaG,BetaOutput);
        
        Dw = [doutdwG doutdthetaG];
        Dx = doutdalphaG(:,1:N);
    end

    ANN.numn = numnF;
    ANN.f = f;
    ANN.df = df;
    ANN.ddf = ddf;
    ANN.BetaOutput = BetaOutput;        
    modelclass.eval_sensitivity_dfdaPen = @eval_sensitivity_dfdaPen;
    function [Da,Daw] = eval_sensitivity_dfdaPen(paramsF,x,a)
        if constrained
            error('not implemented')
        end
        ANN.w = paramsF(1:nwFw);
        ANN.theta = paramsF(nwFw+1:nwFw+nwFtheta);
        [doutdinput,ddoutdinputdw,ddoutdinputdtheta] = ANNeval(ANN,[reshape(x,[],1); a],{'i','iw','it'},'F');
        Da = permute(doutdinput(:,nU+N+1:nU+N+N_alpha),[3,4,1,2]);
        Daw = permute(cat(3, ddoutdinputdw(:,nU+N+1:nU+N+N_alpha,:), ddoutdinputdtheta(:,nU+N+1:nU+N+N_alpha,:)),[4,5,1,2,3]);    
%                 [~,alphaF,betaF] = EvaluateANN(numnF,wF,thetaF,[reshape(dfda_pt(iS,iPt,:),[],1); Alpha(:,iS)],f,BetaOutput);
%                 [doutdalphaF,~,~] = ANNSensitivityBP(numnF,wF,df,alphaF,betaF,BetaOutput);
%                 ResDfda(iS,iPt,:,:) = doutdalphaF(:,nU+N+1:nU+N*N_alpha);         
    end

    modelclass.eval_sensitivity_origStabPen = @eval_sensitivity_origStabPen;
    function [outputF,Dx,Dw,Dwx] = eval_sensitivity_origStabPen(paramsF,x,u)
        if constrained
            error('not implemented')
        end
        [outputF,alphaF,betaF] = EvaluateANN(numnF,paramsF(1:nwFw),paramsF(nwFw+1:nwFw+nwFtheta),[u;x],f,BetaOutput);
        [doutdinput,doutdw,doutdtheta,ddoutdinputdw,ddoutdinputdtheta] = ANNSensitivityFP(numnF,paramsF(1:nwFw),df,ddf,alphaF,betaF,BetaOutput);
        Dx  = doutdinput(:,nU+1:end);
        Dw  = [doutdw doutdtheta];
        Dwx = cat(3,ddoutdinputdw,ddoutdinputdtheta);
    end

    modelclass.eval_sensitivity_steady_state_incr = @eval_sensitivity_steady_state_incr;
    function [Du,Dx,Duw,Dxw] = eval_sensitivity_steady_state_incr(paramsF,x,u,a)
        if constrained
            error('not implemented')
        end
        ANN.w = paramsF(1:nwFw);
        ANN.theta = paramsF(nwFw+1:nwFw+nwFtheta);
        [doutdinput,ddoutdinputdw,ddoutdinputdtheta] = ANNeval(ANN,[u;x;a],{'i','iw','it'},'F');
        idxs_u = 1:nU;
        idxs_x = nU+1:nU+N;
        Du = permute(doutdinput(:,idxs_u),[3,4,1,2]);
        Dx = permute(doutdinput(:,idxs_x),[3,4,1,2]);
        Duw = permute(cat(3, ddoutdinputdw(:,idxs_u,:), ddoutdinputdtheta(:,idxs_u,:)),[3,1,2,4,5]);  
        Dxw = permute(cat(3, ddoutdinputdw(:,idxs_x,:), ddoutdinputdtheta(:,idxs_x,:)),[3,1,2,4,5]);    
    end

    if ~constrained
        modelclass.get_smart_handle_f = @get_smart_handle_f;
        modelclass.get_smart_handle_g = @get_smart_handle_g;
    end
    function fhandle = get_smart_handle_f(paramsF,alpha)
        if nargin < 2
            alpha = [];
        end
        wF = paramsF(1:nwFw);
        thetaF = paramsF(nwFw+1:nwFw+nwFtheta);
        numw = numnF(1:end-1).*numnF(2:end);
        numwcum = [0;cumsum(numw)];
        numncum = cumsum(numnF);
        if length(numnF)==5
            W1 = reshape(wF(numwcum(1)+1:numwcum(2)),numnF(2),numnF(1));
            W2 = reshape(wF(numwcum(2)+1:numwcum(3)),numnF(3),numnF(2));
            W3 = reshape(wF(numwcum(3)+1:numwcum(4)),numnF(4),numnF(3));
            W4 = reshape(wF(numwcum(4)+1:numwcum(5)),numnF(5),numnF(4));
            T1 = thetaF((numncum(1)+1:numncum(2))-numnF(1));
            T2 = thetaF((numncum(2)+1:numncum(3))-numnF(1));
            T3 = thetaF((numncum(3)+1:numncum(4))-numnF(1));
            T4 = thetaF((numncum(4)+1:numncum(5))-numnF(1));
            if problem.nU > 0
                fhandle = @(x,u) W4*f(W3*f(W2*f(W1*[u;x;alpha]-T1)-T2)-T3)-T4;
            else
                fhandle = @(x,u) W4*f(W3*f(W2*f(W1*[x;alpha]-T1)-T2)-T3)-T4;
            end
        elseif length(numnF)==4
            W1 = reshape(wF(numwcum(1)+1:numwcum(2)),numnF(2),numnF(1));
            W2 = reshape(wF(numwcum(2)+1:numwcum(3)),numnF(3),numnF(2));
            W3 = reshape(wF(numwcum(3)+1:numwcum(4)),numnF(4),numnF(3));
            T1 = thetaF((numncum(1)+1:numncum(2))-numnF(1));
            T2 = thetaF((numncum(2)+1:numncum(3))-numnF(1));
            T3 = thetaF((numncum(3)+1:numncum(4))-numnF(1));
            if problem.nU > 0
                fhandle = @(x,u) W3*f(W2*f(W1*[u;x;alpha]-T1)-T2)-T3;
            else
                fhandle = @(x,u) W3*f(W2*f(W1*[x;alpha]-T1)-T2)-T3;
            end
        elseif length(numnF)==3
            W1 = reshape(wF(numwcum(1)+1:numwcum(2)),numnF(2),numnF(1));
            W2 = reshape(wF(numwcum(2)+1:numwcum(3)),numnF(3),numnF(2));
            T1 = thetaF((numncum(1)+1:numncum(2))-numnF(1));
            T2 = thetaF((numncum(2)+1:numncum(3))-numnF(1));
            if problem.nU > 0
                fhandle = @(x,u) W2*f(W1*[u;x;alpha]-T1)-T2;
            else
                fhandle = @(x,u) W2*f(W1*[x;alpha]-T1)-T2;
            end
        else
            error('number of layers not implemented')
        end
    end

    function ghandle = get_smart_handle_g(paramsG)
        wG = paramsG(1:nwGw);
        thetaG = paramsG(nwGw+1:nwGw+nwGtheta);
        numw = numnG(1:end-1).*numnG(2:end);
        numwcum = [0;cumsum(numw)];
        numncum = cumsum(numnG);
        if length(numnG)==4
            W1 = reshape(wG(numwcum(1)+1:numwcum(2)),numnG(2),numnG(1));
            W2 = reshape(wG(numwcum(2)+1:numwcum(3)),numnG(3),numnG(2));
            W3 = reshape(wG(numwcum(3)+1:numwcum(4)),numnG(4),numnG(3));
            T1 = thetaG((numncum(1)+1:numncum(2))-numnG(1));
            T2 = thetaG((numncum(2)+1:numncum(3))-numnG(1));
            T3 = thetaG((numncum(3)+1:numncum(4))-numnG(1));
            ghandle = @(x) W3*f(W2*f(W1*x-T1)-T2)-T3;
        elseif length(numnG)==3
            W1 = reshape(wG(numwcum(1)+1:numwcum(2)),numnG(2),numnG(1));
            W2 = reshape(wG(numwcum(2)+1:numwcum(3)),numnG(3),numnG(2));
            T1 = thetaG((numncum(1)+1:numncum(2))-numnG(1));
            T2 = thetaG((numncum(2)+1:numncum(3))-numnG(1));
            ghandle = @(x) W2*f(W1*x-T1)-T2;
        else
            error('number of layers not implemented')
        end
    end

    if constrained_alpha_single
       modelclass.get_alpha_to_alpha_handler = @get_alpha_to_alpha_handler_alphacon1;
    elseif constrained_alpha_double
       modelclass.get_alpha_to_alpha_handler = @get_alpha_to_alpha_handler_alphacon2;
    end
    
    function hdl = get_alpha_to_alpha_handler_alphacon1(HFmod)
        hdl = @(a) HFmod.f_alpha(x_star_orig,u_star_orig,a);
    end
    function hdl = get_alpha_to_alpha_handler_alphacon2(HFmod)
        hdl = @(a) [HFmod.f_alpha(x_star_orig,u_star_orig,a); HFmod.f_alpha(x_star_orig_2,u_star_orig_2,a)];
    end


    if constrained_alpha_single
       modelclass.fix_alpha = @fix_alpha_alphacon1;
    elseif constrained_alpha_double
       modelclass.fix_alpha = @fix_alpha_alphacon2;
%     TODO: implement?
%     elseif constrained_steadystate
%        modelclass.fix_alpha = @fix_alpha_steadystate;
    end
    function Alpha_0 = fix_alpha_alphacon1(Alpha_0)
        if ~isempty(initalphadir)
            ANNmod = read_model_fromfile(problem,initalphadir);
            for iS = 1:size(Alpha_0,2)
                Alpha_0(:,iS) = ANNmod.f_alpha(x_star_orig,u_star_orig,ANNmod.alpha_learned(iS))./alpha_norm;
            end
        end
    end
    function Alpha_0 = fix_alpha_alphacon2(Alpha_0)
        if ~isempty(initalphadir)
            ANNmod = read_model_fromfile(problem,initalphadir);
            for iS = 1:size(Alpha_0,2)
                Alpha_0(:,iS) = [ANNmod.f_alpha(x_star_orig  ,u_star_orig  ,ANNmod.alpha_learned(iS)) ; ...
                                 ANNmod.f_alpha(x_star_orig_2,u_star_orig_2,ANNmod.alpha_learned(iS))]./alpha_norm;
            end
        end
    end	
%     function Alpha_0 = fix_alpha_steadystate(Alpha_0)
%         %TODO
%     end	

    modelclass.visualize = @visualize;
    function visualize(paramsF,paramsG,opt_vis)
        
        opt_vis.dummy = 0;
        if ~isfield(opt_vis,'label_in')
            opt_vis.label_in = [];
        end
        if ~isfield(opt_vis,'label_out')
            opt_vis.label_out = [];
        end
            
        if constrained
            paramsF = [paramsF;zeros(N,1)];
        end
        wF_curr = paramsF(1:nwFw);
        thetaF_curr = paramsF(nwFw+1:nwFw+nwFtheta);
        
        figure()

        opt_vis.plot_labels_IO = 1;
        if isempty(opt_vis.label_in)
            for i = 1:nU
                opt_vis.label_in{i} = sprintf('{u}_{%d}',i);
            end
            for i = 1:numnF(1)-nU
                iPoint = nU+i;
                opt_vis.label_in{iPoint} = sprintf('{x}_{%d}',i);
            end
        end

        if isempty(opt_vis.label_out)
            for i = 1:numnF(end)
                opt_vis.label_out{i} = sprintf('\\dot{x}_{%d}',i);
            end
        end

        ANN_visualize(numnF,wF_curr,thetaF_curr,opt_vis)
    end

    modelclass.save_compact_file = @save_compact_file;
    function save_compact_file(path,paramsF,paramsG,normalize,x0,opt)
        
        opt.dummy = 0;
        if ~isfield(opt,'mat_file')
            opt.mat_file = 1;
        end
%         if constrained_ref
%             paramsF = [paramsF;zeros(N,1)];
%         end
%         
%         wF = paramsF(1:nwFw);
%         thetaF = paramsF(nwFw+1:nwFw+nwFtheta);
%         
%         if constrained_ref
%             thetaF(end-N+1:end) = EvaluateANN(numnF,wF,thetaF,[u_star;x_star],f,BetaOutput);
%         end
%         
%         A_u = (normalize.u_max-normalize.u_min)/2;
%         B_u = (normalize.u_max+normalize.u_min)/2;
%         A_y = (normalize.y_max-normalize.y_min)/2;
%         B_y = (normalize.y_max+normalize.y_min)/2;
%         for iU = 1:nU
%             [wF,thetaF] = ANNmod_affinetransf_input(numnF,wF,thetaF,A_u(iU),B_u(iU),iU);
%         end
%         for iY = 1:nY
%             [wF,thetaF] = ANNmod_affinetransf_input(numnF,wF,thetaF,A_y(iY),B_y(iY),nU + iY);
%             [wF,thetaF] = ANNmod_affinetransf_output(numnF,wF,thetaF,A_y(iY),0,iY);
%         end

        raw_mod = get_raw_model(paramsF,paramsG,normalize);
        wF = raw_mod.wF;
        thetaF = raw_mod.thetaF;

        numwF = numnF(1:end-1).*numnF(2:end);
        numwcumF = [0;cumsum(numwF)];
        numncumF = cumsum(numnF);
        nLayF = length(numnF);
        N = numnF(end);
        for i = 1:nLayF-1
            W{i} = reshape(wF(numwcumF(i)+1:numwcumF(i+1)),numnF(i+1),numnF(i));
            T{i} = thetaF((numncumF(i)+1:numncumF(i+1))-numnF(1));
        end
        
        if useG
            numwG = numnG(1:end-1).*numnG(2:end);
            numwcumG = [0;cumsum(numwG)];
            numncumG = cumsum(numnG);
            nLayG = length(numnG);
            for i = 1:nLayG-1
                W_G{i} = reshape(wG(numwcumG(i)+1:numwcumG(i+1)),numnG(i+1),numnG(i));
                T_G{i} = thetaG((numncumG(i)+1:numncumG(i+1))-numnG(1));
            end
        end
        
        if opt.mat_file
            save(path,'N', 'nU', 'nY', 'W','T','x0','useG')
            if useG
                save(path,'W_G','T_G','-append')
            end
        else
            create_directory_if_not_found(path);
            % writing setup file
            f_setup = fopen([path '/setup'],'w');
            fprintf(f_setup, '%% number of variables:\n');
            fprintf(f_setup, '%d\n', N);
            fprintf(f_setup, '%% number of hidden layers:\n');
            fprintf(f_setup, '%d\n', nLayF - 2);
            fclose(f_setup);
            
            % writing initial state
            writematrix(x0', [path '/initial_state.csv']);
            
            for i = 1:nLayF-1
                writematrix(W{i} , sprintf('%s/weights_%d.csv',path,i-1));
                writematrix(T{i}', sprintf('%s/biases_%d.csv',path,i-1));
            end
        end
        fprintf('compact file saved!\n')
    end

    modelclass.get_raw_model = @get_raw_model;
    function raw_mod = get_raw_model(paramsF,paramsG,normalize)
        
        if constrained_ref
            paramsF = [paramsF;zeros(N,1)];
        end
        
        wF = paramsF(1:nwFw);
        thetaF = paramsF(nwFw+1:nwFw+nwFtheta);
        
        if constrained_ref
            thetaF(end-N+1:end) = EvaluateANN(numnF,wF,thetaF,[u_star;x_star],f,BetaOutput);
        end
        
        A_u = (normalize.u_max-normalize.u_min)/2;
        B_u = (normalize.u_max+normalize.u_min)/2;
        A_y = (normalize.y_max-normalize.y_min)/2;
        B_y = (normalize.y_max+normalize.y_min)/2;
        for iU = 1:nU
            [wF,thetaF] = ANNmod_affinetransf_input(numnF,wF,thetaF,A_u(iU),B_u(iU),iU);
        end
        if ~useG
            for iY = 1:nY
                [wF,thetaF] = ANNmod_affinetransf_input(numnF,wF,thetaF,A_y(iY),B_y(iY),nU + iY);
                [wF,thetaF] = ANNmod_affinetransf_output(numnF,wF,thetaF,A_y(iY),0,iY);
            end
        end
        
        for iX = 1:N
            [wF,thetaF] = ANNmod_affinetransf_output(numnF,wF,thetaF,1/normalize.t_norm,0,iX);
        end
        
        if useG
            wG = paramsG(1:nwGw);
            thetaG = paramsG(nwGw+1:nwGw+nwGtheta);
            for iY = 1:nY
                [wG,thetaG] = ANNmod_affinetransf_output(numnG,wG,thetaG,A_y(iY),B_y(iY),iY);
            end
        end 

        raw_mod.useG = useG;
        raw_mod.numnF = numnF;
        raw_mod.wF = wF;
        raw_mod.thetaF = thetaF;
        if useG
            raw_mod.numnG = numnG;
            raw_mod.wG = wG;
            raw_mod.thetaG = thetaG;
        end
    end

end