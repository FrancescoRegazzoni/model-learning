function mod = read_model_fromfile(problem,model_folder,opt)
    % Loads a pre-trained model from the folder 'model_folder'.

    %% options
    opt.dummy = 0;
    if ~isfield(opt,'idxNet')
        opt.idxNet = Inf;
    end
    if ~isfield(opt,'show_history')
        opt.show_history = 0;
    end
    if ~isfield(opt,'get_history')
        opt.get_history = 0;
    end

    %% File paths
    directory = problem.dir_nets;
    disp('=======================================================')
    fprintf('%s/%s --- iter %d\n',directory,model_folder,opt.idxNet);
    disp('=======================================================')
    fprintf('loading ... ');

    baseopt = get_base_options();
    fullANNpath = strcat(baseopt.BaseDir,'/',directory,'/',model_folder);
    mod.full_path = fullANNpath;
    optionsfile = strcat(fullANNpath,'/options.ini');
    fileOpt = strcat(fullANNpath,'/data.mat');
    NormalizFile = strcat(fullANNpath,'/normalization.mat');
    if isinf(opt.idxNet)
        fileNet = strcat(fullANNpath,'/net.mat');
    else
        fileNet = strcat(fullANNpath,'/networks/',sprintf('net_%06d.mat',opt.idxNet));
    end

    %% Options from file
    N                   = iniread(optionsfile,'Model','N');
    N_alpha             = iniread(optionsfile,'Model','N_alpha','d',0);
    useG                = iniread(optionsfile,'Model','useG','d',0);
    modelclassname      = iniread(optionsfile,'Model','modelclass','s','@modelclass_ANN');
    dt_integration      = iniread(optionsfile,'Numerics','dt_integration','s');
    
    
    ds_def.type         = iniread(optionsfile,'Problem','dataset_type','s','file');
    ds_def.source_train = iniread(optionsfile,'Problem','dataset_source_train','s');
    ds_def.source_tests = iniread(optionsfile,'Problem','dataset_source_tests','s');
    ds_def.source       = iniread(optionsfile,'Problem','dataset_source','s');
    ds_def.problem = problem;

    netw = load(fileNet);
    opts = load(fileOpt,'xhat_init');

    if isfield(netw,'paramsF')
        paramsF = netw.paramsF;
    else 
        % kept for legacy versions
        paramsF = [netw.wF; netw.thetaF];
        opts2 = load(fileOpt,'N','N_alpha');
        N = opts2.N;
        if isfield(opts2, 'N_alpha')
            N_alpha = opts2.N_alpha;
        else
            N_alpha = 0;
        end
    end
    paramsG = [];
    if exist(NormalizFile,'file')
        normFile = load(NormalizFile);
        normalize = normFile.normalization;
        if isfield(normalize,'alpha_to_alpha')
            mod.alpha_to_alpha = normalize.alpha_to_alpha;
        end
    else
        % kept for legacy versions
        normalize = load(fileOpt,'u_min','u_max','y_min','y_max','t_norm');
        if ~isfield(normalize,'u_min' ), normalize.u_min  = -1; end
        if ~isfield(normalize,'u_max' ), normalize.u_max  =  1; end
        if ~isfield(normalize,'y_min' ), normalize.y_min  = -1; end
        if ~isfield(normalize,'y_max' ), normalize.y_max  =  1; end
        if ~isfield(normalize,'t_norm'), normalize.t_norm =  1; end

        normalize = adapt_dimension_struct_field(normalize,'u_min',problem.nU);
        normalize = adapt_dimension_struct_field(normalize,'u_max',problem.nU);
        normalize = adapt_dimension_struct_field(normalize,'y_min',problem.nY);
        normalize = adapt_dimension_struct_field(normalize,'y_max',problem.nY);

    end
    mod.paramsF = paramsF;

    modeclasshandler = eval(modelclassname);
    modelclass = modeclasshandler(optionsfile,problem,N,N_alpha,useG);

    %% renormalization  


    A_u = (normalize.u_max-normalize.u_min)/2;
    B_u = (normalize.u_max+normalize.u_min)/2;
    A_y = (normalize.y_max-normalize.y_min)/2;
    B_y = (normalize.y_max+normalize.y_min)/2;
    if ~useG
        A_x = [A_y; ones(N-problem.nY,1)];
        B_x = [B_y; zeros(N-problem.nY,1)];
    else
        A_x = ones(N,1);
        B_x = zeros(N,1);
    end
    if ~isfield(normalize,'alpha_norm')
        normalize.alpha_norm = ones(N_alpha,1);
    end
    if isstruct(normalize.alpha_norm)
        A_a = (normalize.alpha_norm.max-normalize.alpha_norm.min)/2;
        B_a = (normalize.alpha_norm.max+normalize.alpha_norm.min)/2;
    else
        A_a = normalize.alpha_norm;
        B_a = zeros(N_alpha,1);
    end   
    

    f_renormalized = 0;
    g_renormalized = 0;

    % u_normalization = @(u) (u-B_u)./A_u;
    % y_normalization = @(y) (y-B_y)./A_y;
    % t_normalization = @(t) t / t_norm;
    if isfield(modelclass,'apply_normalization')
    %     normaliz.u_norm = u_normalization;
    %     normaliz.y_norm = y_normalization;
    %     normaliz.t_norm = t_normalization;
        modelclass.apply_normalization(normalize);
    end

    x0 = opts.xhat_init;
    if ~useG
        x0(1:problem.nY) = A_y.*x0(1:problem.nY) + B_y;
    end

    if isfield(modelclass,'set_xref')
        modelclass.set_xref(opts.xhat_init);
    end

    %% Model definition
    mod.problem = problem;
    mod.nX = N;
    % mod.N_alpha = N_alpha;
    mod.nA = N_alpha;

    mod.x0 = x0; 
    if dt_integration(1) == '@'
        error('not implemented')
    end
    mod.dt = str2double(dt_integration);   

    mod.advance_type = 'nonlinear_explicit';
    if N_alpha > 0
    %     mod.particularize_model = @particularize_model;
        f_alpha_base = @(x,u,a) modelclass.eval_f(x,u,a,paramsF,mod.dt);
        if f_renormalized
            mod.f_alpha = f_alpha_base;
        else
            mod.f_alpha = @(x,u,a) A_x.*f_alpha_base((x-B_x)./A_x,(u-B_u)./A_u,(a-B_a)./A_a)/normalize.t_norm;
        end
    else
        if isfield(modelclass,'get_smart_handle_f')
            mod.f = get_f_normalized(modelclass.get_smart_handle_f(paramsF));
        else
            mod.f = get_f_normalized(@(x,u) modelclass.eval_f(x,u,[],paramsF,mod.dt));
        end
    end
    % function mod_part = particularize_model(mod,alpha)
    %     mod_part = mod;
    %     mod_part.alpha = alpha;
    %     if isfield(modelclass,'get_smart_handle_f')
    %         mod_part.f = get_f_normalized(modelclass.get_smart_handle_f(paramsF,alpha));
    %     else
    %         mod_part.f = get_f_normalized(@(x,u) modelclass.eval_f(x,u,alpha,paramsF,mod.dt));
    %     end
    % end

    function f = get_f_normalized(f_base)
        if f_renormalized
            f = f_base;
        else
            f = @(x,u) A_x.*f_base((x-B_x)./A_x,(u-B_u)./A_u)/normalize.t_norm;
        end
    end

    if problem.samples_variability 
        samples = load(fileOpt,'d');
        if isfield(samples.d.train{1}, 'alpha')
            alpha_original = [];
            for iS = 1:length(samples.d.train)
               alpha_original = [alpha_original samples.d.train{iS}.alpha];
            end
            mod.alpha_original = alpha_original;
        end
    end

    if N_alpha > 0
        mod.alpha_learned = B_a + netw.Alpha.*A_a;
        mod.alpha_min = min(mod.alpha_learned,[],2);
        mod.alpha_max = max(mod.alpha_learned,[],2);
        for iA = 1:N_alpha
            if mod.alpha_min(iA) == mod.alpha_max(iA)
                if mod.alpha_min(iA) == 0
                    mod.alpha_min = -1;
                    mod.alpha_max = 1;
                else
                    mod.alpha_min = mod.alpha_min * 0.9;
                    mod.alpha_max = mod.alpha_max * 1.1;
                end
            end        
        end
    end
    % if isfield(netw,'IC')
    if ~problem.fixed_x0
        mod.IC_learned = A_x.*reshape(netw.IC,N,[]) + B_x;
        mod.x0_min = min(mod.IC_learned,[],2);
        mod.x0_max = max(mod.IC_learned,[],2);
        for iX = 1:N
            if mod.x0_min(iX) == mod.x0_max(iX)
                if mod.x0_min(iX) == 0
                    mod.x0_min = -1;
                    mod.x0_max = 1;
                else
                    mod.x0_min = mod.x0_min * 0.9;
                    mod.x0_max = mod.x0_max * 1.1;
                end
            end        
        end
    end 
    mod.u_implicit = 0;

    if useG
        mod.output_type = 'nonlinear';
        paramsG = netw.paramsG;
        mod.paramsG = netw.paramsG;
        if isfield(modelclass,'get_smart_handle_g')
            g_base = modelclass.get_smart_handle_g(mod.paramsG);
        else
            g_base = @(x) modelclass.eval_g(x,netw.paramsG);
        end
        if g_renormalized
            mod.g = g_base;
        else
            mod.g = @(x) A_y.*g_base((x-B_x)./A_x) + B_y;
        end
    else
        mod.output_type = 'insidestate'; 
    end

    dfdx_base = @dfdx;
    if f_renormalized
        mod.dfdx = dfdx_base;
    else
        mod.dfdx = @(x,u,a) A_x.*dfdx_base((x-B_x)./A_x,(u-B_u)./A_u,(a-B_a)./A_a)./A_x'/normalize.t_norm;
    end
    function Dx = dfdx(x,u,a)
        modelclass.eval_f(x,u,a,paramsF,mod.dt);
        [~,Dx,~] = modelclass.eval_sensitivity_f();
    end

    dfda_base = @dfda;
    if f_renormalized
        mod.dfda = dfda_base;
    else
        mod.dfda = @(x,u,a) A_x.*dfda_base((x-B_x)./A_x,(u-B_u)./A_u,(a-B_a)./A_a)./A_a'/normalize.t_norm;
    end
    function Da = dfda(x,u,a)
        modelclass.eval_f(x,u,a,paramsF,mod.dt);
        [~,~,Da] = modelclass.eval_sensitivity_f();
    end

    if useG
        dgdx_base = @dgdx;
        if g_renormalized
            mod.dgdx = dgdx_base;
        else
            mod.dgdx = @(x) A_y.*dgdx_base((x-B_x)./A_x)./A_x';
        end
    end
    function Dx = dgdx(x)  
        modelclass.eval_g(x,paramsG); 
        [~,Dx] = eval_sensitivity_g() ;
    end

    if isfield(modelclass,'get_alpha_to_alpha_handler')
        mod.get_alpha_to_alpha_handler = modelclass.get_alpha_to_alpha_handler;    
    end

    mod.NetworkDir = model_folder;
    mod.idxNet = opt.idxNet;

    if opt.get_history
        mod.EE_n = netw.EE_n;
        mod.EE_n_test = netw.EE_n_test;
    end
    
    mod.datasets_def = ds_def;

    if isfield(modelclass,'visualize')
        mod.visualize = @visualize;    
    end
    function visualize(opt_vis)
        opt_vis.dummy = 0;
        modelclass.visualize(paramsF,paramsG,opt_vis);
    end

    if isfield(modelclass,'save_compact_file')
        mod.save_compact_file = @save_compact_file;
        mod.save_compact_dir = @save_compact_dir;
    end
    function save_compact_file()
        modelclass.save_compact_file([fullANNpath,'/compact.mat'],paramsF,paramsG,normalize,x0,struct('mat_file',1));
    end
    function save_compact_dir()
        modelclass.save_compact_file([fullANNpath,'/compact'],paramsF,paramsG,normalize,x0,struct('mat_file',0));
    end

    if isfield(modelclass,'get_raw_model')
        mod.get_raw_model = @() modelclass.get_raw_model(paramsF,paramsG,normalize);
    end


    if opt.show_history
        figure()

        if ~isfield(netw,'Act')
            netw.Act = ones(1,size(netw.EE_w,1));
        end

        subplot(1,2,1)
        loglog(sqrt(netw.EE_t)','k--','linewidth',2);
        hold on
        loglog(sqrt(netw.EE_w(netw.Act==1,:))','linewidth',2);
        ax = gca;
        ax.ColorOrderIndex = 1;
        loglog(sqrt(netw.EE_w_test)','-','linewidth',1);
        grid on              
        if isfield(netw,'Leg')
    %         legend([{'ERR'} netw.Leg{netw.Act==1} {'EFFdiff (test)'}] , ...
    %                'location','southwest');
            format_legend = @(s,v) sprintf('%s (%1.2e)',s,v);
            legend_labels{1} = format_legend('ERR',sqrt(netw.EE_t(end)));
            for iE = 1:length(netw.Act)
                if netw.Act(iE)
                    legend_labels = [legend_labels {format_legend(netw.Leg{iE},sqrt(netw.EE_w(iE,end)))}];
                end
            end
            legend_labels = [legend_labels {format_legend('EFFdiff (test)',sqrt(netw.EE_w_test(end)))}];
            legend(legend_labels,'location','southwest');
        end
        title('weighted values')

        subplot(1,2,2)
        loglog(sqrt(netw.EE_n(netw.Act==1,:))','linewidth',2);
        hold on
        ax = gca;
        ax.ColorOrderIndex = 1;
        loglog(sqrt(netw.EE_n_test)','-','linewidth',1);
        grid on
        title('absolute values')
    end

    fprintf('done!\n');

end