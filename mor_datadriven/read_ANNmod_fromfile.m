function ANNmod = read_ANNmod_fromfile(problem,model_folder,opt)
    % Loads a pre-trained model from the folder 'model_folder'.
    % This is a legacy version of read_model_fromfile. This function should be
    % used only to load models trained with old versions of the library.

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

directory = problem.dir_nets;

disp('=======================================================')
fprintf('%s/%s --- iter %d\n',directory,model_folder,opt.idxNet);
disp('=======================================================')
fprintf('loading ... ');

baseopt = get_base_options();
fullANNpath = strcat(baseopt.BaseDir,'/',directory,'/',model_folder);
fileNetwork = strcat(fullANNpath,'/network.mat');
fileOpt = strcat(fullANNpath,'/data.mat');
if isinf(opt.idxNet)
    fileNet = strcat(fullANNpath,'/net.mat');
else
    fileNet = strcat(fullANNpath,'/networks/',sprintf('net_%06d.mat',opt.idxNet));
end
if exist(fileNetwork,'file')
    % Old-style output
	netw = load(fileNetwork);
    
    numnF = netw.numnF;
    N = netw.N;
    useG = netw.useG;
    f = netw.f;
    if useG
        numnG = netw.numnG;
        wG = netw.wG;
        thetaG = netw.thetaG; 
    end
    
    wF = netw.wF;
    thetaF = netw.thetaF;
else
    % New-style output
	netw = load(fileNet);
%     warning('off','MATLAB:load:variableNotFound');
    warning off
	opts = load(fileOpt,'numnF','N','N_alpha','useG','f','df','xhat_init','dt_integration',...
        'u_min','u_max','y_min','y_max','t_norm','int_dyn');
    warning on
%     warning('on','MATLAB:load:variableNotFound');

    numnF = opts.numnF;
    N = opts.N;
    useG = opts.useG;
    if isfield(opts, 'N_alpha')
        N_alpha = opts.N_alpha;
    else
        N_alpha = 0;
    end
    if isfield(opts, 'int_dyn')
        int_dyn = opts.int_dyn;
    else
        int_dyn = 1;
    end
    f = opts.f;
    df = opts.df;
    if useG
        optsG = load(fileOpt,'numnG');
        numnG = optsG.numnG;
        wG = netw.wG;
        thetaG = netw.thetaG; 
    end
    if problem.samples_variability
        samples = load(fileOpt,'d');
        alpha_original = [];
        for iS = 1:length(samples.d.train)
           alpha_original = [alpha_original samples.d.train{iS}.alpha];
        end
    end
    
    if ~int_dyn
        mod_ord = N/problem.nY;
    end
    
    wF = netw.wF;
    thetaF = netw.thetaF;  
    
    x0 = opts.xhat_init;
        
    if opts.dt_integration(1) == '@'
        error('not implemented')
    end
    dt = str2double(opts.dt_integration);

    if ~isfield(opts,'u_min' ), opts.u_min  = -1; end
    if ~isfield(opts,'u_max' ), opts.u_max  =  1; end
    if ~isfield(opts,'y_min' ), opts.y_min  = -1; end
    if ~isfield(opts,'y_max' ), opts.y_max  =  1; end
    if ~isfield(opts,'t_norm'), opts.t_norm =  1; end
    
    % renormalization    
    if length(opts.u_min)==1
        opts.u_min = opts.u_min*ones(problem.nU,1);
    end
    if length(opts.u_max)==1
        opts.u_max = opts.u_max*ones(problem.nU,1);
    end
    if length(opts.y_min)==1
        opts.y_min = opts.y_min*ones(problem.nY,1);
    end
    if length(opts.y_max)==1
        opts.y_max = opts.y_max*ones(problem.nY,1);
    end    
    for iU = 1:problem.nU
        A = (opts.u_max(iU)-opts.u_min(iU))/2;
        B = (opts.u_max(iU)+opts.u_min(iU))/2;
        [wF,thetaF] = ANNmod_affinetransf_input(numnF,wF,thetaF,A,B,iU);
    end
    for iY = 1:problem.nY
        A = (opts.y_max(iY)-opts.y_min(iY))/2;
        B = (opts.y_max(iY)+opts.y_min(iY))/2;
        if ~useG
            [wF,thetaF] = ANNmod_affinetransf_input(numnF,wF,thetaF,A,B,problem.nU + iY);
            if ~int_dyn
                for iOrd = 2:mod_ord
                    [wF,thetaF] = ANNmod_affinetransf_input(numnF,wF,thetaF,A,B,problem.nU + iY + (iOrd-1)*problem.nY);
                end
            end
            [wF,thetaF] = ANNmod_affinetransf_output(numnF,wF,thetaF,A,0,iY);
            x0(iY) = A*x0(iY) + B;
%             x0(iY) = x0(iY) * (opts.y_max(iY)-opts.y_min(iY))/2 + (opts.y_max(iY)+opts.y_min(iY))/2;
        end
        if useG
            [wG,thetaG] = ANNmod_affinetransf_output(numnG,wG,thetaG,A,B,iY);
        end
    end
    for iOut = 1:numnF(end)
        [wF,thetaF] = ANNmod_affinetransf_output(numnF,wF,thetaF,1/opts.t_norm,0,iOut);
    end
    
    if opt.get_history
        ANNmod.EE_n = netw.EE_n;
        ANNmod.EE_n_test = netw.EE_n_test;
    end

    if opt.show_history
        figure()

        subplot(1,2,1)
        loglog(sqrt(netw.EE_t)','k--','linewidth',2);
        hold on
        loglog(sqrt(netw.EE_w)','linewidth',2);
        ax = gca;
        ax.ColorOrderIndex = 1;
        loglog(sqrt(netw.EE_w_test)','--','linewidth',2);
        % legend('ERR'            , ...
               % 'ERRdiff'        , ...
               % 'ERRratio'       , ...
               % 'ERRorigEq'      , ...
               % 'ERRorigStab'    , ...
               % 'ERRf'           , ...
               % 'ERRg'           , ...
               % 'ERRend'         , ...
               % 'ERRraise'       , ...
               % 'EFFdiff (test)' , ...
               % 'location','southwest');
        grid on                

        subplot(1,2,2)
        loglog(sqrt(netw.EE_n)','linewidth',2);
        hold on
        ax = gca;
        ax.ColorOrderIndex = 1;
        loglog(sqrt(netw.EE_n_test)','--','linewidth',2);
        grid on
        pause(1e-16)
    end
end

%% smart evaluation of f
numw = numnF(1:end-1).*numnF(2:end);
numwcum = [0;cumsum(numw)];
numncum = cumsum(numnF);
if length(numnF)==4

    W1 = reshape(wF(numwcum(1)+1:numwcum(2)),numnF(2),numnF(1));
    W2 = reshape(wF(numwcum(2)+1:numwcum(3)),numnF(3),numnF(2));
    W3 = reshape(wF(numwcum(3)+1:numwcum(4)),numnF(4),numnF(3));
    T1 = thetaF((numncum(1)+1:numncum(2))-numnF(1));
    T2 = thetaF((numncum(2)+1:numncum(3))-numnF(1));
    T3 = thetaF((numncum(3)+1:numncum(4))-numnF(1));
    %T3(1) = T3(1) + .002;

    fANN = @(x) W3*f(W2*f(W1*x-T1)-T2)-T3;
    if problem.nU > 0
        ANNmod.f = @(x,u) W3*f(W2*f(W1*[u;x]-T1)-T2)-T3;
    else
        ANNmod.f = @(x,u) W3*f(W2*f(W1*x-T1)-T2)-T3;
    end
elseif length(numnF)==3

    W1 = reshape(wF(numwcum(1)+1:numwcum(2)),numnF(2),numnF(1));
    W2 = reshape(wF(numwcum(2)+1:numwcum(3)),numnF(3),numnF(2));
    T1 = thetaF((numncum(1)+1:numncum(2))-numnF(1));
    T2 = thetaF((numncum(2)+1:numncum(3))-numnF(1));

    fANN = @(x) W2*f(W1*x-T1)-T2;
    if problem.nU > 0
        ANNmod.f = @(x,u) W2*f(W1*[u;x]-T1)-T2;
    else
        ANNmod.f = @(x,u) W2*f(W1*x-T1)-T2;
    end
else
    error('number of layers not implemented')
end

%% smart evaluation of g
if useG
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
        %T3(1) = T3(1) + .002;

        gANN = @(x) W3*f(W2*f(W1*x-T1)-T2)-T3;
    elseif length(numnG)==3

        W1 = reshape(wG(numwcum(1)+1:numwcum(2)),numnG(2),numnG(1));
        W2 = reshape(wG(numwcum(2)+1:numwcum(3)),numnG(3),numnG(2));
        T1 = thetaG((numncum(1)+1:numncum(2))-numnG(1));
        T2 = thetaG((numncum(2)+1:numncum(3))-numnG(1));

        gANN = @(x) W2*f(W1*x-T1)-T2;
    else
        error('number of layers not implemented')
    end
end

%% Model definition

ANNmod.numnF   = numnF  ;
ANNmod.wF      = wF     ;
ANNmod.thetaF  = thetaF ;
ANNmod.fANN = fANN;
ANNmod.int_dyn = int_dyn;

ANNmod.useG    = useG   ;
if useG
    ANNmod.numnG   = numnG  ;
    ANNmod.wG      = wG     ;
    ANNmod.thetaG  = thetaG ;
    ANNmod.gANN = gANN;
end
ANNmod.N_alpha = N_alpha;
if problem.samples_variability
    ANNmod.alpha_original = alpha_original;
end
if N_alpha > 0
    ANNmod.alpha_learned = netw.Alpha;
end
if isfield(netw,'IC')
    ANNmod.IC_learned = reshape(netw.IC,N,[]);
end
ANNmod.problem = problem;
ANNmod.nX = N;
ANNmod.x0 = x0;
ANNmod.dt = dt;
ANNmod.u_implicit = 0;

ANNmod.advance_type = 'nonlinear_explicit';
% if problem.nU > 0
%     ANNmod.f = @(x,u) ANN([u;x]);
% else
%     ANNmod.f = @(x,u) ANN(x);
% end
if useG
    ANNmod.output_type = 'nonlinear'; 
    ANNmod.g = gANN;
else
    ANNmod.output_type = 'insidestate'; 
end

ANNmod.dfdx = @dfdx;
ANNmod.dfdalpha = @dfdalpha;

% ANNmod.advance_type = 'ANN';
% ANNmod.output_type = 'insidestate'; 

% ANNmod.advance_type = 'f';
% ANNmod.f = @(x,u) ANN([u;x]);
% % ANNmod.f = @(x,u) W2*f(W1*[u;x]-T1)-T2;
% % ANNmod.output_type = 'linear'; 
% % ANNmod.G = [eye(problem.nY) zeros(problem.nY,N-problem.nY)];
% % ANNmod.g0 = zeros(problem.nY,1);
% ANNmod.output_type = 'nonlinear';
% nY = problem.nY;
% ANNmod.g = @(x) x(1:nY);

ANNmod.NetworkDir = model_folder;
ANNmod.idxNet = opt.idxNet;

fprintf('done!\n');

function diff = dfdx(x)        
    [~,alphaF,betaF] = EvaluateANN(numnF,wF,thetaF,x,f,1);
    [doutdalphaF,~,~] = ANNSensitivityBP(numnF,wF,df,alphaF,betaF,1);
    diff = doutdalphaF(:,problem.nU+1:problem.nU+N);
end
function diff = dfdalpha(x)        
    [~,alphaF,betaF] = EvaluateANN(numnF,wF,thetaF,x,f,1);
    [doutdalphaF,~,~] = ANNSensitivityBP(numnF,wF,df,alphaF,betaF,1);
    diff = doutdalphaF(:,problem.nU+N+1:problem.nU+N+N_alpha);
end

end