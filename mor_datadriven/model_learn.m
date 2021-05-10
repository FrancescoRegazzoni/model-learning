function model_learn(optionsfile)

% Learns a model form a collection of input-output pairs. All the options
% are passed by an option file (see mor_datadriven/opt_example.ini for a
% full list of available options and their documentation).

%% Fixed Options
checkGradOption = 0; % true se si vuole controllare che i gradienti calcolati via LM siano coerenti
xNormPlot = 1;
xmean_liminf = .1;
xmean_limsup = 2;
xmean_target = .4;
normalizeplot = 0;

%% Load options
example_name        = iniread(optionsfile,'Problem','example','s');
problem_file        = iniread(optionsfile,'Problem','problem','s');
ds_def.type         = iniread(optionsfile,'Problem','dataset_type','s','file');
ds_def.source_train = iniread(optionsfile,'Problem','dataset_source_train','s');
ds_def.source_tests = iniread(optionsfile,'Problem','dataset_source_tests','s');
ds_def.source       = iniread(optionsfile,'Problem','dataset_source','s');
noise_y             = iniread(optionsfile,'Problem','noise_y','d',0);
bias_y              = iniread(optionsfile,'Problem','bias_y','d',0);
noise_u             = iniread(optionsfile,'Problem','noise_u','d',0);
bias_u              = iniread(optionsfile,'Problem','bias_u','d',0);

N                   = iniread(optionsfile,'Model','N');
N_alpha             = iniread(optionsfile,'Model','N_alpha','d',0);
useG                = iniread(optionsfile,'Model','useG','d',0);
xhat_init           = iniread(optionsfile,'Model','x_init','d',0)*ones(N,1);
xhat_end            = iniread(optionsfile,'Model','x_end','d',0)*ones(N,1);
modelclassname      = iniread(optionsfile,'Model','modelclass','s','@modelclass_ANN');
    
%DirOut              = iniread(optionsfile,'Output','dir','s');
Name                = iniread(optionsfile,'Output','name','s');
saveNetwork         = iniread(optionsfile,'Output','save','d',1);
add_time_stamp      = iniread(optionsfile,'Output','add_time_stamp','d',1);
saveHistory         = iniread(optionsfile,'Output','history','d',1);
savePeriod          = iniread(optionsfile,'Output','period_save','s','1');
snapPeriod          = iniread(optionsfile,'Output','period_snap','s','@(k) (k<10) || ((k < 500) && mod(k,10)==0) || mod(k,100)==0');
plotPeriod          = iniread(optionsfile,'Output','period_plot','s','@(k) (k<10) || ((k < 500) && mod(k,10)==0) || mod(k,100)==0');

% u_min               = iniread(optionsfile,'Normalization','u_min','d',-1);
% u_max               = iniread(optionsfile,'Normalization','u_max','d',1);
% y_min               = iniread(optionsfile,'Normalization','y_min','d',-1);
% y_max               = iniread(optionsfile,'Normalization','y_max','d',1);
u_min               = iniread(optionsfile,'Normalization','u_min','d',-inf);
u_max               = iniread(optionsfile,'Normalization','u_max','d',inf);
y_min               = iniread(optionsfile,'Normalization','y_min','d',-inf);
y_max               = iniread(optionsfile,'Normalization','y_max','d',inf);
t_norm              = iniread(optionsfile,'Normalization','t_norm','d',1);

x0_min              = iniread(optionsfile,'Initialization','x0_min','d',0);
x0_max              = iniread(optionsfile,'Initialization','x0_max','d',0);
alpha_min           = iniread(optionsfile,'Initialization','alpha_min','d',0);
alpha_max           = iniread(optionsfile,'Initialization','alpha_max','d',0);

customPenalizations = iniread(optionsfile,'Penalizations','custom_penalizations','s','');
diffPen_s           = iniread(optionsfile,'Penalizations','pen_diff','s','1');
diffPen_type        = iniread(optionsfile,'Penalizations','pen_diff_type','s','q');
diffPen_discr_toll  = iniread(optionsfile,'Penalizations','pen_diff_discrete_tol','d',.5);
diffPen_discr_coeff = iniread(optionsfile,'Penalizations','pen_diff_discrete_coeff_within_int','d',1e-2);
alphaPen_s          = iniread(optionsfile,'Penalizations','pen_end','s','0');
alphaRaise_s        = iniread(optionsfile,'Penalizations','pen_raise','s','0');
targetRaise         = iniread(optionsfile,'Penalizations','RaiseTarget');
ratioPen_s          = iniread(optionsfile,'Penalizations','pen_ratio','s','0');
ratio_type          = iniread(optionsfile,'Penalizations','ratio_type','d',2);
ratio_pen_out       = iniread(optionsfile,'Penalizations','ratio_penalize_outputs','d',0);

algorithm           = iniread(optionsfile,'Numerics','algorithm','s','LM');
LinSearchAlg        = iniread(optionsfile,'Numerics','linearsearch','d',2);
nmax_iter           = iniread(optionsfile,'Numerics','nmax_iter','d',inf);
ComputeGradientsNum = iniread(optionsfile,'Numerics','ComputeGradientsNum','d',0);
epNum               = iniread(optionsfile,'Numerics','epNum','d',1e-4);
SA                  = iniread(optionsfile,'Numerics','SA','d',0);
SA_stepsize_s       = iniread(optionsfile,'Numerics','SA_stepsize','s');
SA_T_s              = iniread(optionsfile,'Numerics','SA_T','s');
%perEval             = iniread(optionsfile,'Numerics','perEval');
dt_integration      = iniread(optionsfile,'Numerics','dt_integration','s');
dt_evaluation       = iniread(optionsfile,'Numerics','dt_evaluation','s','x1');
% interpolation_mode_u= iniread(optionsfile,'Numerics','interpolation_mode_u','s','mean_forward');
interpolation_mode_u= iniread(optionsfile,'Numerics','interpolation_mode_u','s','pointwise');
interpolation_mode_y= iniread(optionsfile,'Numerics','interpolation_mode_y','s','pointwise');

do_plot             = iniread(optionsfile,'plots','do_plot','d',1);
plot_showU          = iniread(optionsfile,'plots','show_u','d',1);

problem = problem_get(example_name,problem_file);
baseopt = get_base_options();
BaseDir = baseopt.BaseDir;
SaveDir = strcat(BaseDir,'/',problem.dir_nets,'/');

exampledir = [baseopt.ExampleDir '/' example_name];
cd(exampledir);
       
ds_def.problem = problem;
d = datasetcouple_get(ds_def);

% modeclasshandler = eval(['@' modelclassname]);
modeclasshandler = eval(modelclassname);
modelclass = modeclasshandler(optionsfile,problem,N,N_alpha,useG);
 
%% Check
if problem.samples_variability && N_alpha == 0
    warning('Attention! The problem allows for samples variability, but N_alpha = 0')
end
if ~problem.samples_variability && N_alpha > 1
    warning('Attention! N_alpha > 0, but the problem does not allow for samples variability')
end

if useG
    ratio_pen_out = 1; % otherwise the first nY variables would be excluded from ratio penalization
end

%% Datasets preprocessing (time interpolation, noise)
nS = length(d.train);
if isfield(d,'tests')
    nS_test = length(d.tests);
else
    nS_test = 0;
end
% nU = size(d.train{1}.uu,1);
% nY = size(d.train{1}.yy,1);
nU = problem.nU;
nY = problem.nY;

noise_y = adapt_dimension(noise_y,nY);
bias_y = adapt_dimension(bias_y,nY);
noise_u = adapt_dimension(noise_u,nU);
bias_u = adapt_dimension(bias_u,nU);
is_y_noised_man = 0;
if sum(abs(noise_y)) + sum(abs(bias_y)) > 0
    is_y_noised_man = 1;
end
is_u_noised_man = 0;
if sum(abs(noise_u)) + sum(abs(bias_u)) > 0
    is_u_noised_man = 1;
end

is_y_noised = 0;
if isfield(d.train{1},'yy_ex') % assuming tht if the first sample is noised, al the database is so (and viceversa)
    is_y_noised = 1;
end
if is_y_noised && is_y_noised_man
    warning('Adding noise to a yet noised database');
end
is_y_noised = is_y_noised || is_y_noised_man;

for jS = 1:nS
       
    if is_y_noised_man
        d.train{jS}.yy_ex = d.train{jS}.yy;
    end
    if nU > 0 && is_u_noised_man
        d.train{jS}.uu_ex = d.train{jS}.uu;
    end
    
    tr{jS} = mor_ANN_interpolate_test(problem,d.train{jS},dt_integration,dt_evaluation,interpolation_mode_u,interpolation_mode_y);
    
    % noise + bias
    rng('default') % to ensure reproducibility
    ddt = [tr{jS}.dt(1) .5*(tr{jS}.dt(2:end)+tr{jS}.dt(1:end-1)) tr{jS}.dt(end)];
    if is_y_noised_man
        tr{jS}.yy = tr{jS}.yy + noise_y./sqrt(ddt).*randn(size(tr{jS}.yy)) + bias_y;
    end
    if nU > 0 && is_u_noised_man
        tr{jS}.uu = tr{jS}.uu + noise_u./sqrt(ddt).*randn(size(tr{jS}.uu)) + bias_u;
    end
end

for jS = 1:nS_test
    ts{jS} = mor_ANN_interpolate_test(problem,d.tests{jS},dt_integration,dt_evaluation,interpolation_mode_u,interpolation_mode_y);
end

num_train_closed_loop = 0;
for jS = 1:nS
    if tr{jS}.closed_loop
        num_train_closed_loop = num_train_closed_loop + 1;
    end
end

%% Renormalization

u_min = adapt_dimension(u_min,nU);
u_max = adapt_dimension(u_max,nU);
y_min = adapt_dimension(y_min,nY);
y_max = adapt_dimension(y_max,nY);

for iU = 1:nU
    if u_min(iU) == -inf
        u_min(iU) = problem.u_min(iU);
    end
    if u_max(iU) == inf
        u_max(iU) = problem.u_max(iU);
    end
end
for iY = 1:nY
    if y_min(iY) == -inf
        y_min(iY) = problem.y_min(iY);
    end
    if y_max(iY) == inf
        y_max(iY) = problem.y_max(iY);
    end
end

normalization.u_min = u_min;
normalization.u_max = u_max;
normalization.y_min = y_min;
normalization.y_max = y_max;
normalization.t_norm = t_norm;
%     u_normalization = @(u) (2*u - u_min - u_max)./(u_max-u_min);
u_normalization = @(u) renormalize_m1p1(u,u_min,u_max);
% y_normalization = @(y) (2*y - y_min - y_max)./(y_max-y_min);
y_normalization = @(y) renormalize_m1p1(y,y_min,y_max);
t_normalization = @(t) t / t_norm;

if isfield(modelclass,'apply_normalization')
    normalization.alpha_norm = modelclass.apply_normalization(normalization);
end
alpha_to_alpha_handler = [];
if isfield(modelclass,'get_alpha_to_alpha_handler')
    alpha_to_alpha_handler = modelclass.get_alpha_to_alpha_handler;
end

for jS = 1:nS
    tr{jS} = NormalizeTest(tr{jS});
end
for jS = 1:nS_test
    ts{jS} = NormalizeTest(ts{jS});
end
if problem.fixed_x0
    if isfield(problem,'y_ref')
        y_ref = y_normalization(problem.y_ref);
    else
        % if y_ref is not available, it is estimated as the mean value of
        % the inital values of the trainset
        y_ref = zeros(nY,1);
        for jS = 1:nS
            y_ref = y_ref + tr{jS}.yy(:,1);
        end
        y_ref = y_ref/nS;
    end
end

function out = NormalizeTest(in)
    out = in;
    out.tt = t_normalization(in.tt);
    out.dt = t_normalization(in.dt);
    out.dtEval = t_normalization(in.dtEval);
    out.Tmax = t_normalization(in.Tmax);
    if nU > 0
        out.uu = u_normalization(in.uu);
    end
    out.yy = y_normalization(in.yy);
    if isfield(in,'yy_ex')
        out.yy_ex = y_normalization(in.yy_ex);
    end
end

%% General Initialization

if ~useG
    if problem.fixed_x0
        xhat_init(1:nY) = y_ref;
        xhat_end(1:nY)  = y_ref;
    end
end

    diffPen_act = ~( strcmp(    diffPen_s,'0') || isempty(    diffPen_s));   
   alphaPen_act = ~( strcmp(   alphaPen_s,'0') || isempty(   alphaPen_s));
 alphaRaise_act = ~( strcmp( alphaRaise_s,'0') || isempty( alphaRaise_s));
   ratioPen_act = ~( strcmp(   ratioPen_s,'0') || isempty(   ratioPen_s));

    diffPen_fun = inline(    diffPen_s,'k');
   alphaPen_fun = inline(   alphaPen_s,'k');
 alphaRaise_fun = inline( alphaRaise_s,'k');
   ratioPen_fun = inline(   ratioPen_s,'k');

% if (size(uu,2)~=nT || size(uu_test,2)~=nT || size(yy,2)~=nT || size(yy_test,2)~=nT ...
%         || size(yy,3) ~= nS || size(yy_test,3) ~= nS_test)
% if (size(yy,3) ~= nS || size(yy_test,3) ~= nS_test)    
%     error('Non consistent dimensions')
% end

nTevalTot = 0;
for iS = 1:nS
    nTevalTot = nTevalTot+tr{iS}.nTeval;
end

normL2tot_train = 0;
normL2tot_test = 0;
% y_zero_normalized = y_normalization(zeros(nY,1));
y_zero_normalized = y_normalization(problem.y0_norm_computation);
for iS = 1:nS
    if tr{iS}.eval_diff
        normL2tot_train = normL2tot_train+sum(sum((tr{iS}.yy(:,tr{iS}.idxEval) - y_zero_normalized).^2,1).*tr{iS}.dtEval);
    end
end
for iS = 1:nS_test
    normL2tot_test = normL2tot_test+sum(sum((ts{iS}.yy(:,ts{iS}.idxEval) - y_zero_normalized).^2,1).*ts{iS}.dtEval);
end
   
nwF = modelclass.nwF;
if useG
    nwG = modelclass.nwG;
else
    nwG = 0;
end
ntotdof = nwF + nwG;

savePeriod_fun = get_period_function(savePeriod);
snapPeriod_fun = get_period_function(snapPeriod);
plotPeriod_fun = get_period_function(plotPeriod);

function fun = get_period_function(period_string)
    if period_string(1) == '@'
        fun = eval(period_string);
    else
        fun = @(k) mod(k,str2num(period_string)) == 0;
    end
end

if isfield(modelclass,'set_xref')
    modelclass.set_xref(xhat_end);
end

switch diffPen_type
    case 'q'
        loss_diff = @(x) x;
        d_loss_diff = @(x) ones(length(x),1);
    case 'd'
        toll_normalized = diffPen_discr_toll ./ ((normalization.y_max - normalization.y_min)/2);
        loss_diff = @(x) diffPen_discr_coeff*min(toll_normalized,max(-toll_normalized,x)) ...
                         + max(0,x - toll_normalized) ...
                         + min(0,x + toll_normalized);
        d_loss_diff = @(x) diffPen_discr_coeff + (1-diffPen_discr_coeff)*(abs(x) >= toll_normalized);
    otherwise
        error(['Penalization type ' diffPen_type ' not recognized'])
end


%% Initial conditions initialization

if ~problem.fixed_x0
    nIC = N*nS;
    x0_min = adapt_dimension(x0_min,N);
    x0_max = adapt_dimension(x0_max,N);
%         if ~useG
%             if ~isfield(problem,'y0_min'), problem.y0_min = problem.y_min; end
%             if ~isfield(problem,'y0_max'), problem.y0_max = problem.y_max; end
%             x0_min(1:nY) = y_normalization(problem.y0_min);
%             x0_max(1:nY) = y_normalization(problem.y0_max);
%         end
    IC_0 = x0_min + (x0_max-x0_min).*rand(N,nS);
    if ~useG
        % if not useG, the initial guess for the first nY initial conditions are the measured outputs
        for iS = 1:nS, IC_0(1:nY,iS) = tr{iS}.yy(:,1); end        
    end
    % if necessary, correct the initial conditions
    if isfield(modelclass,'fix_IC')
        IC_0 = modelclass.fix_IC(IC_0,tr);
    end
else
    nIC = 0;
end

%% Samples-specific parameters initialization (alpha)

alpha_manually_fixed = 0;
if N_alpha > 0
    alpha_min = adapt_dimension(alpha_min,N_alpha);
    alpha_max = adapt_dimension(alpha_max,N_alpha);
    Alpha_0 = alpha_min + (alpha_max-alpha_min).*rand(N_alpha,nS);
    nAL = N_alpha*nS;
    if isfield(modelclass,'fix_alpha')
        Alpha_0 = modelclass.fix_alpha(Alpha_0);
        alpha_manually_fixed = 1;
    end
else
    nAL = 0;
end

%% penalizations
pen_handlers{1} = penalization_params_direct(optionsfile);
pen_handlers{2} = penalization_origin_equilibrium(optionsfile);
pen_handlers{3} = penalization_dfdalpha(optionsfile);
pen_handlers{4} = penalization_f_constrained_wrt_alpha(optionsfile);
pen_handlers{5} = penalization_ref_equilibrium(optionsfile);
pen_handlers{6} = penalization_alpha_equilibria(optionsfile);
pen_handlers{7} = penalization_steady_state_increasing(optionsfile);
pen_handlers{8} = penalization_rhs_increasing(optionsfile);
pen_handlers{9} = penalization_initial_state_equilibrium(optionsfile);

custom_pen_list = strsplit(customPenalizations,';');
for ip = 1:length(custom_pen_list)
    if ~isempty(custom_pen_list{ip})
        penalization_handler = eval(custom_pen_list{ip});
        pen_handlers = [pen_handlers {penalization_handler(optionsfile)}];
    end
end

misc.xhat_end = xhat_end;
misc.normaliz = normalization;
misc.tr = tr;
for ipen = 1:length(pen_handlers)
    for ipenloc=1:pen_handlers{ipen}.num_pen
        coef_string = pen_handlers{ipen}.coefficient_string{ipenloc};
        pen_handlers{ipen}.pen_active(ipenloc) = ~( strcmp(coef_string,'0') || isempty(coef_string)); 
        pen_handlers{ipen}.pen_fun{ipenloc}.handler = inline(coef_string,'k');
    end
    pen_handlers{ipen} = pen_handlers{ipen}.initialize(pen_handlers{ipen},problem,N,N_alpha,useG,nS,nwF,nwG,nAL,nIC,misc);
    if isfield(pen_handlers{ipen},'alpha_norm')
        if isfield(normalization,'alpha_norm')
            if ~isempty(normalization.alpha_norm) && ~isequal(normalization.alpha_norm, pen_handlers{ipen}.alpha_norm)
                error('different alpha_norm specified by multiple agents (modelclass and/or penalizations')
            end
        end
        normalization.alpha_norm = pen_handlers{ipen}.alpha_norm;
    end
    
    if isfield(pen_handlers{ipen},'get_alpha_to_alpha_handler')
        if isempty(alpha_to_alpha_handler)
            alpha_to_alpha_handler = pen_handlers{ipen}.get_alpha_to_alpha_handler;
        else
            error('alpha_to_alpha_handler specified by multiple agents (modelclass and/or penalizations')
        end
    end

    if N_alpha > 0
        if any(pen_handlers{ipen}.pen_active) && isfield(pen_handlers{ipen},'fix_alpha')
            Alpha_0_new = pen_handlers{ipen}.fix_alpha(pen_handlers{ipen},problem,Alpha_0);
            if alpha_manually_fixed
                if isequal(Alpha_0_new,Alpha_0)
                    error('Alpha_0 set to different values by multiple agents')
                end
            end
            Alpha_0 = Alpha_0_new;
            alpha_manually_fixed = 1;
        end
    end
end

alpha_norm_abstent = 0;
if ~isfield(normalization,'alpha_norm')
    alpha_norm_abstent = 1;
elseif isempty(normalization.alpha_norm) 
    alpha_norm_abstent = 1;
end
if alpha_norm_abstent
    normalization.alpha_norm = ones(N_alpha,1);
end

if ~isempty(alpha_to_alpha_handler)
    HFmod = problem.get_model(problem);
    normalization.alpha_to_alpha = alpha_to_alpha_handler(HFmod);
end

%% Oputput initialization

maxRows = 100;
ratioOutputTestTrain=.3;
ratioHistory=.3;
marginHistoryLR = .02;
marginHistoryTB = .05;
marginPlot = .001;

EE_t = []; % total
EE_n = []; % non weighted
EE_w = []; % weighted
EE_n_test = [];
EE_w_test = [];
if do_plot
    FigMain = figure('units','normalized','outerposition',[0 0 1 1]);
end
numIter = 0;
% stichastic gradient descent
idx_SamplesSelection = [];

if nU > 0
    uMin = inf*ones(nU,1);
    uMax = -inf*ones(nU,1);
    for jS = 1:nS
        uMin = min(uMin,min(tr{jS}.uu,[],2));
        uMax = max(uMax,max(tr{jS}.uu,[],2));
    end
    for jS = 1:nS_test
        uMin = min(uMin,min(ts{jS}.uu,[],2));
        uMax = max(uMax,max(ts{jS}.uu,[],2));
    end
end

if do_plot
    yBounds = [-1.1 1.1];

    nStotal = nS+nS_test;

    if nStotal <= maxRows^2
        nRows = ceil(sqrt(nStotal));
        nSOut_test = nS_test;
        nSOut = nS;
    else
        nRows = maxRows;
        nSOut_test = round(ratioOutputTestTrain*nRows^2);
        nSOut = nRows^2-nSOut_test;
        if (nSOut_test > nS_test)
            nSOut_test = nS_test;
            nSOut = nRows^2-nSOut_test;
        elseif(nSOut > nS)
            nSOut = nS;
            nSOut_test = nRows^2-nSOut;
        end
    end
    nCols = nRows;    
    axHistory_1 = subplot('Position',[marginHistoryLR .5+marginHistoryTB ratioHistory-2*marginHistoryLR .5-2*marginHistoryTB]);
    axHistory_2 = subplot('Position',[marginHistoryLR marginHistoryTB ratioHistory-2*marginHistoryLR .5-2*marginHistoryTB]);

    widthPlot = (1-ratioHistory)/nCols;
    heigthPlot = 1/nRows;
    colCurr = 0;
    rowCurr = 0;
    for j = 1:nSOut
        leftCurr = ratioHistory + rowCurr*widthPlot;
        bottomCurr = 1 - (colCurr+1)*heigthPlot;
        ax_train(j) = subplot('Position',[leftCurr+marginPlot bottomCurr+marginPlot widthPlot-2*marginPlot heigthPlot-2*marginPlot]);

        rowCurr = rowCurr+1;
        if rowCurr == nRows
            rowCurr = 0;
            colCurr = colCurr+1;
        end
    end

    for j = 1:nSOut_test
        leftCurr = ratioHistory + rowCurr*widthPlot;
        bottomCurr = 1 - (colCurr+1)*heigthPlot;
        ax_test(j) = subplot('Position',[leftCurr+marginPlot bottomCurr+marginPlot widthPlot-2*marginPlot heigthPlot-2*marginPlot]);

        rowCurr = rowCurr+1;
        if rowCurr == nRows
            rowCurr = 0;
            colCurr = colCurr+1;
        end
    end

end

%% Saving initialization
if saveNetwork
    % base directory
    if ~exist(SaveDir,'dir')
        if ~mkdir(SaveDir)
            error('unable to create directory')
        end
    end
    CurrentDirName = strcat(Name,modelclass.name_modelclass,'_N',num2str(N));
    if isfield(modelclass,'laystrF')
        CurrentDirName = [CurrentDirName modelclass.laystrF];
    end
    if useG && isfield(modelclass,'laystrG')
        CurrentDirName = [CurrentDirName modelclass.laystrG];
    end
    CurrentDirName = [CurrentDirName '_dof' num2str(ntotdof)];
    CurrentDirName = [CurrentDirName '_ntrain' num2str(nS)];
    if add_time_stamp
        CurrentDirName = [CurrentDirName '_' datestr(now,'yyyy-mm-dd_hh-MM-ss')];
    end
    CurrentDir = strcat(SaveDir,CurrentDirName);
    mkdir(SaveDir,CurrentDirName);
    fprintf('created directory %s\n',CurrentDirName);
    try
        clipboard('copy',CurrentDirName)
    catch
    end
    % copy files
    copyfile(strcat(mfilename('fullpath'),'.m'),CurrentDir)
    copyfile(optionsfile,CurrentDir)
    copyfile(optionsfile,[CurrentDir '/options.ini'])
    copyfile(problem_file,CurrentDir)
    % output files
    OutputFile = strcat(CurrentDir,'/data.mat');
    DiaryFile = strcat(CurrentDir,'/diary.txt');
    NormalizFile = strcat(CurrentDir,'/normalization.mat');
    diary(DiaryFile)
    save(OutputFile,'-regexp','^(?!(FigMain)$).') % exclude FigMain
    save(NormalizFile,'normalization') % exclude FigMain
    if isfield(modelclass,'save_to_file_init'), modelclass.save_to_file_init(OutputFile); end    
    fileNet = strcat(CurrentDir,'/net.mat');
    
    if saveHistory
        mkdir(CurrentDir,'networks');
    end
end

%% Printing out machine info
[~, hostname] = system('hostname');
hostname = strtrim(hostname);
pid = feature('getpid');

fprintf('Time: %s\n',datestr(now,'YYYY/mm/dd HH:MM:SS.FFF'));
fprintf('Host: %s\n',hostname);
fprintf('PID:  %d\n\n',pid);

%% Optimization loop

clear d

numIter_last = -1;

if ~isfield(modelclass,'paramsF0')
    modelclass.paramsF0 = randn(modelclass.nwF,1);
end
x0 = modelclass.paramsF0;
if useG
    if ~isfield(modelclass,'paramsG0')
        modelclass.paramsG0 = randn(modelclass.nwG,1);
    end
    x0 = [x0;modelclass.paramsG0];
end
if ~problem.fixed_x0
    x0 = [x0;reshape(IC_0,nIC,1)];
end
if N_alpha > 0
    x0 = [x0;reshape(Alpha_0,nAL,1)];
end

switch upper(algorithm)
    case 'SD'
        alphaFix = .1;
        SteepestDescent(x0,@FuncGrad,LinSearchAlg,@PlotFnc,alphaFix,nmax_iter);
    case 'LBFGS'
        lbfgs(x0,@Func,@Grad,@PlotFnc,nmax_iter);
    case 'LM'
        opt.nMaxIter = nmax_iter;
        opt.LSopt.algorithm = LinSearchAlg;
        opt.LSopt.alphaFix = .1;
        %opt.dphiZeroThreshold = 1e-5;
        opt.dphiZeroThreshold = -1; %negative to deactivate
        opt.SA = SA;
        if SA
            opt.SA_stepsize = inline(SA_stepsize_s,'k');
            opt.SA_T = inline(SA_T_s,'k');
        end
        if ComputeGradientsNum
            LevenbergMarquardt(x0,@Func,@FuncGradNum,@LSFuncGradNum,@PlotFnc,opt);
        else
            LevenbergMarquardt(x0,@Func,@FuncGrad,@LSFuncGrad,@PlotFnc,opt);
        end
    case 'TEST'
        TestDerivativesNumerically(x0,@Evaluation);
    case 'TEST_F'
        modelclass_test_derivatives(modelclass,N,N_alpha,nU);
    otherwise
        error('Algorithm not recognized')
end

fprintf('finished!\n')

%% Functions definition

function ret = Func(x)
    ret = Evaluation(x,1,0,0,0,0);
end

function ret = Grad(x)
    ret = Evaluation(x,0,1,0,0,0);
end

function xNew = PlotFnc(x)
    xNew = Evaluation(x,0,0,0,0,1);
end

function [E,DE] = FuncGrad(x)
    if nargout == 1
        E = Evaluation(x,1,0,0,0,0);
    else
        [E,DE] = Evaluation(x,1,1,0,0,0);
    end
end

function [F,DF] = LSFuncGrad(x)
    [F,DF] = Evaluation(x,0,0,1,1,0);
end

function [E,DE] = FuncGradNum(x)
    E = Evaluation(x,1,0,0,0,0);
    nx = nwF+nwG;
    DE = zeros(nx,1);
    for ix = 1:nx
        xnew = x;
        xnew(ix) = xnew(ix) + epNum;
        Enew = Evaluation(xnew,1,0,0,0,0);
        DE(ix) = (Enew-E)/epNum;
    end
end

function [F,DF] = LSFuncGradNum(x)
    F = Evaluation(x,0,0,1,0,0);
    nx = nwF+nwG;
    nF = length(F);
    DF = zeros(nF,nx);
    for ix = 1:nx
        xnew = x;
        xnew(ix) = xnew(ix) + epNum;
        Fnew = Evaluation(xnew,0,0,1,0,0);
        DF(:,ix) = (Fnew-F)/epNum;
    end
end

function [varargout] = Evaluation(Params,compErr,compGrad,compRes,compJac,step_finalization)
           
    checkGrad = checkGradOption;
    if checkGrad && compJac
        compGrad = 1;
    else
        checkGrad = 0;
    end
     
    compGradRaise = 1;
    if ~(compGrad || checkGrad) || ~alphaRaise_act
        compGradRaise = 0;
    end
    
    if step_finalization
        fprintf('**************************************************\n')
        if saveNetwork
            fprintf('Output directory: %s\n',CurrentDirName)
        end
    end
    
    idxCurr = 0;
    paramsF =      Params(idxCurr+1:idxCurr+nwF); idxCurr = idxCurr+nwF;
    if useG
        paramsG =  Params(idxCurr+1:idxCurr+nwG); idxCurr = idxCurr+nwG;
    else
        paramsG = [];
    end
    if ~problem.fixed_x0
        IC =       Params(idxCurr+1:idxCurr+nIC); idxCurr = idxCurr+nIC;
        IC = reshape(IC,N,nS);
    else
        IC = [];
    end
    if N_alpha > 0
        Alpha =    Params(idxCurr+1:idxCurr+nAL); idxCurr = idxCurr+nAL;
        Alpha = reshape(Alpha,N_alpha,nS);
    else
        Alpha = [];
    end
    
    for ipen = 1:length(pen_handlers)
        for ipenloc=1:pen_handlers{ipen}.num_pen
            pen_handlers{ipen}.coef(ipenloc) = pen_handlers{ipen}.pen_fun{ipenloc}.handler(numIter);
        end
    end

        diffPen  =     diffPen_fun(numIter);
       alphaPen  =    alphaPen_fun(numIter);
     alphaRaise  =  alphaRaise_fun(numIter);
       ratioPen  =    ratioPen_fun(numIter);
    
    if num_train_closed_loop == 0
        ratioPen = 0;
    end
    
    if numIter > numIter_last
        % it means that this is the first evaluation for this optimization epoch
        
%         idx_SamplesSelection = sort(randperm(nS,nS));
        idx_SamplesSelection = sort(randperm(nS,round(nS)));
        
        for ipen = 1:length(pen_handlers)
            if isfield(pen_handlers{ipen},'init_epoch')
                pen_handlers{ipen}.init_epoch(pen_handlers{ipen},paramsF,paramsG,Alpha,IC);
            end
        end
        
        numIter_last = numIter;
    end

    if step_finalization
        nStot = nS+nS_test;
        %uutot = cat(3,uu,uu_test);
    else
        nStot = nS;
        %uutot = uu;
    end
    
    xmean_mean = zeros(N,1);
    xmeanL2_mean = zeros(N,1);
    
    if ratio_pen_out
        nVarRatio = N;
        idxVarRatio = 1:N;
    else
        nVarRatio = N-nY;
        idxVarRatio = nY+1:N;
    end
    
%     normFactDiff = 1/(nS*nY*tt(end));
%     normFactDiff_test = 1/(nS_test*nY*tt(end));    
%     normFactDiff = 1/sum(sum((sum((yy(:,idxEval,1:nS)).^2,1).*dtEval)));
%     normFactDiff_test = 1/sum(sum((sum((yy_test(:,idxEval,1:nS_test)).^2,1).*dtEval)));   
    normFactDiff = 1/normL2tot_train; 
    normFactDiff_test = 1/normL2tot_test; 
    normFactEnd = 1/(nS*N);
    normFactRaise = 1/N;
    normFactRatio = 1/(num_train_closed_loop*nVarRatio);
%     normFactRatio = 1/(nS*nVarRatio);
   
    ERRdiff = 0;
    ERRdiff_ex = 0;
    ERRdiff_test = 0;
    ERRend = 0;
    ERRraise = 0;
    ERRratio = 0;
        
    nRes = nY*nTevalTot;
    
    if alphaRaise > 0
        XRaise = zeros(N,1);
    end
    if compGrad
        DG = zeros(nwG,1);
        DF = zeros(nwF,1);
        DI = zeros(nIC,1);
        DA = zeros(nAL,1);
    end
    if compRes
        Res = zeros(nRes,1);
        if alphaPen > 0
            nResEnd = N*nS;
            ResEnd = zeros(nResEnd,1);
        end
%         if alphaRaise > 0
%             ResRaise = -ones(N,1)*sqrt(alphaRaise)*targetRaise;
%         end
        if ratioPen > 0
%             nResRatio = nVarRatio*nS;  
            nResRatio = nVarRatio*num_train_closed_loop;
            ResRatio = zeros(nResRatio,1);
        end
    end
    if compJac    
        JacG = zeros(nRes,nwG);
        JacF = zeros(nRes,nwF);  
        JacI = zeros(nRes,nIC);   
        JacA = zeros(nRes,nAL);   
        if alphaPen > 0
            JacEndG = zeros(nResEnd,nwG); % non va riempita (sono effettivamente zeri)
            JacEndF = zeros(nResEnd,nwF);
            JacEndI = zeros(nResEnd,nIC);
            JacEndA = zeros(nResEnd,nAL);
        end
        if ratioPen > 0
            JacRatioG = zeros(nResRatio,nwG); % non va riempita (sono effettivamente zeri)
            JacRatioF = zeros(nResRatio,nwF);
            JacRatioI = zeros(nResRatio,nIC);
            JacRatioA = zeros(nResRatio,nAL);
        end
    end
    if compJac || compGradRaise
        if alphaRaise > 0
            JacRaiseG = zeros(N,nwG); % non va riempita (sono effettivamente zeri)
            JacRaiseF = zeros(N,nwF);
            JacRaiseI = zeros(N,nIC);
            JacRaiseA = zeros(N,nAL);
        end
    end
     
    if diffPen > 0 || alphaPen > 0 || alphaRaise > 0 || ratioPen > 0
        iResStart = 0;
        iS_closed_loop = 0;
        for iS=1:nStot
            
%             if ismember(iS,idx_SamplesSelection)
%                 skip = 0;
%             else
%                 skip = 1;
%             end
%             if plotta
%                 if plotPeriod_fun(numIter)
%                     skip = 0;
%                 end
%             end
        
            if iS <= nS
                q = tr{iS};
            else
                q = ts{iS - nS};
            end
            if q.closed_loop
                iS_closed_loop = iS_closed_loop+1;
            end
            nTeval = q.nTeval;
            idxEval = q.idxEval;
%             if ~skip
            Tmax = q.Tmax;
            nT = q.nT;
            tt = q.tt;
            dt = q.dt;
            dtEval = q.dtEval;

            %iResStart = (iS-1)*nTeval*nY;

            x = zeros(N,nT);
            z = zeros(N,nT);
            y = zeros(nY,nT);
            e = zeros(nY,nT);
            de = zeros(nY,nT);

            if compGrad || compJac
                Dwf = zeros(N,nwF,nT);
                Dxf = zeros(N,N,nT);
                if N_alpha > 0
                    Daf = zeros(N,N_alpha,nT);
                end
                %Dxg = zeros(nY,N,nT);  
                Dxg = zeros(nY,N,nTeval); 
            end

            % Forward equation
            iTeval = 1;
            for iT = 1:nT
                stepEval = 0;
                if iTeval<=nTeval
                    if idxEval(iTeval) == iT
                        stepEval = 1;
                    end
                end

                if iT == 1
                    if problem.fixed_x0
                        x(:,iT) = xhat_init;
                    else
                        if iS <= nS
                            x(:,iT) = IC(:,iS); 
                        else
                            error('if x0 is not fixed, you cannot have test samples')
                            %x(1:nY,iT) = q.yy(:,1);
                        end
                    end
                else

                    if nU > 0
                        u_curr = q.uu(:,iT-1);
                    else
                        u_curr = [];
                    end
                    if N_alpha > 0
                        a_curr = Alpha(:,iS);
                    else
                        a_curr = [];
                    end
                    x(:,iT) = x(:,iT-1) + dt(iT-1)*modelclass.eval_f(x(:,iT-1),u_curr,a_curr,paramsF,dt(iT-1));

                end

    % 			if iT<nT
                if useG
                    y(:,iT) = modelclass.eval_g(x(:,iT),paramsG);
                else
                    y(:,iT) = x(1:nY,iT);
                end
                if iT<nT
                    if iS<=nS
                        diff = y(:,iT) - q.yy(:,iT);
                        e(:,iT) = loss_diff(diff);
                        if stepEval && q.eval_diff && diffPen > 0
                            if (compGrad || compJac)
                                de(:,iT) = d_loss_diff(diff);
                            end
                            ERRdiff = ERRdiff + dtEval(iTeval)*(e(:,iT)'*e(:,iT));
                            if is_y_noised
                                e_ex = loss_diff(y(:,iT) - q.yy_ex(:,iT));
                                ERRdiff_ex = ERRdiff_ex + dtEval(iTeval)*(e_ex'*e_ex);                                
                            end
                        end
                    else
                        if stepEval
                            e_test = loss_diff(y(:,iT) - q.yy(:,iT));
                            ERRdiff_test = ERRdiff_test + dtEval(iTeval)*(e_test'*e_test);
                        end
                    end
                end
    % 			end

                if iS<=nS
                    if iT>1 && (compGrad || compJac)
                        [Dw,Dx,Da] = modelclass.eval_sensitivity_f();
                        Dwf(:,:,iT-1) = Dw;
                        Dxf(:,:,iT-1) = Dx;
                        if N_alpha > 0
                            Daf(:,:,iT-1) = Da;
                        end
                    end
                    
                    if iT<nT && stepEval
                        if (compGrad || compJac) 
                            if useG
                                [Dw,Dx] = modelclass.eval_sensitivity_g(); 
                                if compGrad 
                                    %TODO: if q.eval_diff?
                                    DG = DG + 2*normFactDiff*dtEval(iTeval)*Dw'*(e(:,iT).*de(:,iT));
                                end
                                Dxg(:,:,iTeval) = Dx;
                            else
                                Dxg(:,:,iTeval) = [eye(nY) zeros(nY,N-nY)];
                            end
                        end

                        if (compRes || compJac) 
                            if q.eval_diff && diffPen > 0
                                iRes =  iResStart + (iTeval-1)*nY;
                                iiRes = iRes+1:iRes+nY;
                                if compJac && useG
                                    JacG(iiRes,:) = sqrt(dtEval(iTeval))*Dw.*de(:,iT);
                                end
                                if compRes
                                    Res(iiRes) = sqrt(dtEval(iTeval))*e(:,iT);
                                end
                            end 
                        end
                    end

                end
                if stepEval
                    iTeval = iTeval+1;
                end
            end

            if iS<=nS

                idxIC = (iS-1)*N+1:iS*N;
                idxAL = (iS-1)*N_alpha+1:iS*N_alpha;
%                 xhat_endCurr = xhat_end;
%                 if ~useG
%                     xhat_endCurr(1:nY) = yy(1:nY,end,iS);
%                 end

                if step_finalization
                    xmean_mean   = xmean_mean   + sum(x(:,1:end-1).*dt,2)/Tmax  - xhat_end;
                    xmeanL2_mean = xmeanL2_mean + sqrt(sum((x(:,1:end-1) - xhat_end).^2.*dt,2)/Tmax);
                end

                if alphaPen > 0
                    Eend = x(:,end) - xhat_end;
                    ERRend = ERRend + (Eend'*Eend);

                    if compRes
                        iResEnd = (iS-1)*N;
                        iiResEnd = iResEnd+1:iResEnd+N;
                        ResEnd(iiResEnd) = Eend;
                    end
                end

                if alphaRaise > 0
                    XRaise = XRaise + sum(x,2);
    %                 if compRes
    %                     ResRaise = ResRaise + sqrt(alphaRaise)*sum(x,2)/nS/nT;
    %                 end                
                end

                if ratioPen > 0 && q.closed_loop
                    switch ratio_type
                        case 1
                            xmean = sum(x(:,1:end-1).*dt,2)/Tmax  - xhat_end;
                        case 2
                            xmean = sqrt(sum((x(:,1:end-1) - xhat_end).^2.*dt,2)/Tmax);
                    end
                    Eend = (x(:,end) - xhat_end);
                    EendNorm = Eend./xmean;                    
                    ERRratio = ERRratio + (EendNorm(idxVarRatio)'*EendNorm(idxVarRatio));

                    if compRes
                        iResEnd = (iS_closed_loop-1)*nVarRatio;
                        iiResEnd = iResEnd+1:iResEnd+nVarRatio;
                        ResRatio(iiResEnd) = EendNorm(idxVarRatio);
                    end    

                end

                %% Gradient
                if compGrad
                    % Backward equation
                    if alphaPen > 0
                        z(:,end) = 2*alphaPen^2*normFactEnd*Eend;
                    else
                        z(:,end) = 0;
                    end
                    if ratioPen > 0 && q.closed_loop
                        zInitRatioPen = 2*ratioPen^2*normFactRatio*EendNorm./xmean;
                        if ~ratio_pen_out
                            zInitRatioPen(1:nY) = 0;
                        end
                        z(:,end) = z(:,end) + zInitRatioPen;
                    end
                    iTeval = nTeval;
                    for iT = nT-1:-1:1
%                         if iT>1 
                            z(:,iT) = z(:,iT+1) + dt(iT) * Dxf(:,:,iT)'*z(:,iT+1);
                            if diffPen > 0 && q.eval_diff
                                if idxEval(iTeval) == iT
                                    %z(:,iT) = z(:,iT) + dtEval(iTeval)*Dxg(:,:,iT)'*e(:,iT);
                                    z(:,iT) = z(:,iT) + 2*diffPen^2*normFactDiff*dtEval(iTeval)*Dxg(:,:,iTeval)'*(e(:,iT).*de(:,iT));
                                    iTeval = iTeval-1;
                                end
                            end
                            if ratioPen > 0 && q.closed_loop
                                switch ratio_type
                                    case 1
                                        zUpdateRatioPen = - 2*ratioPen^2*normFactRatio*Eend.^2.*xmean.^(-3)*dt(iT)/Tmax; 
                                    case 2
                                        %zUpdateRatioPen = - 4*ratioPen^2*normFactRatio*Eend.*xmean.^(-3)*dt(iT)/Tmax.*(x(:,iT) - xhat_end);
                                        zUpdateRatioPen = - 2*ratioPen^2*normFactRatio*Eend.^2.*xmean.^(-4)*dt(iT)/Tmax.*(x(:,iT) - xhat_end);
                                end
                                if ~ratio_pen_out
                                    zUpdateRatioPen(1:nY) = 0;
                                end
                                z(:,iT) = z(:,iT) + zUpdateRatioPen;
                            end
%                         end
                        DF = DF + dt(iT)*Dwf(:,:,iT)'*z(:,iT+1);
                        if N_alpha > 0
                            DA(idxAL) = DA(idxAL) + dt(iT)*Daf(:,:,iT)'*z(:,iT+1);
                        end
                    end
                    if ~problem.fixed_x0
                        DI(idxIC) = z(:,1);
                    end
                end

                %% Jacobians
                if compJac 
                    if diffPen > 0 && q.eval_diff
                        for iY = 1:nY                    
                            for iTeval=1:nTeval                     
                                %iRes = (iS-1)*nTeval*nY + (iTeval-1)*nY + iY; 
                                iRes =  iResStart + (iTeval-1)*nY + iY;                           
                                %z(:,idxEval(iTeval)) = sqrt(dtEval(iTeval))*Dxg(iY,:,idxEval(iTeval));
                                z(:,idxEval(iTeval)) = sqrt(dtEval(iTeval))*Dxg(iY,:,iTeval)*de(iY,idxEval(iTeval));
                                for iT = idxEval(iTeval)-1:-1:1
%                                     if iT>1
                                        z(:,iT) = z(:,iT+1) + dt(iT) * Dxf(:,:,iT)'*z(:,iT+1);
%                                     end
                                    JacF(iRes,:) = JacF(iRes,:) + (dt(iT)*Dwf(:,:,iT)'*z(:,iT+1))';
                                    if N_alpha > 0
                                        JacA(iRes,idxAL) = JacA(iRes,idxAL) + (dt(iT)*Daf(:,:,iT)'*z(:,iT+1))';
                                    end
                                end
                                if ~problem.fixed_x0
                                    JacI(iRes,idxIC) = z(:,1)';
                                end
                            end
                        end
                    end

                    if alphaPen > 0
                        for iX = 1:N 
                            iRes = (iS-1)*N+iX;
                            z(:,end) = zeros(N,1);
                            z(iX,end) = 1;
                            for iT = nT-1:-1:1
%                                 if iT>1
                                    z(:,iT) = z(:,iT+1) + dt(iT) * Dxf(:,:,iT)'*z(:,iT+1);
%                                 end
                                JacEndF(iRes,:) = JacEndF(iRes,:) + (dt(iT)*Dwf(:,:,iT)'*z(:,iT+1))';
                                if N_alpha > 0
                                    JacEndA(iRes,idxAL) = JacEndA(iRes,idxAL) + (dt(iT)*Daf(:,:,iT)'*z(:,iT+1))';
                                end                                
                            end   
                            if ~problem.fixed_x0
                                JacEndI(iRes,idxIC) = z(:,1)'; 
                            end
                        end
                    end

                    if ratioPen > 0 && q.closed_loop
                        for iVar = 1:nVarRatio
                            iX = idxVarRatio(iVar);
                            versorX = zeros(N,1);
                            versorX(iX) = 1;
                            iRes = (iS_closed_loop-1)*nVarRatio+iVar;
                            z(:,end) = zeros(N,1);
                            z(iX,end) = 1/xmean(iX);
                            for iT = nT-1:-1:1
%                                 if iT>1
                                    switch ratio_type
                                        case 1
                                            z(:,iT) = z(:,iT+1) + dt(iT) * Dxf(:,:,iT)'*z(:,iT+1) ...
                                                - xmean(iX)^(-2)*dt(iT)/Tmax*Eend(iX)*versorX;
                                        case 2
                                            z(:,iT) = z(:,iT+1) + dt(iT) * Dxf(:,:,iT)'*z(:,iT+1) ...
                                                - xmean(iX)^(-3)*dt(iT)/Tmax*Eend(iX)*(x(iX,iT) - xhat_end(iX))*versorX;
                                    end                                    
%                                 end
                                JacRatioF(iRes,:) = JacRatioF(iRes,:) + (dt(iT)*Dwf(:,:,iT)'*z(:,iT+1))';
                                if N_alpha > 0
                                    JacRatioA(iRes,idxAL) = JacRatioA(iRes,idxAL) + (dt(iT)*Daf(:,:,iT)'*z(:,iT+1))';
                                end                                
                            end 
                            if ~problem.fixed_x0
                                JacRatioI(iRes,idxIC) = z(:,1)';    
                            end
                        end
                    end
                end
                if compJac || compGradRaise
                    if alphaRaise > 0
                        for iX = 1:N 
                            z(:,end) = zeros(N,1);
                            forcing = zeros(N,1);
                            forcing(iX) = 1/nS/nT;
                            for iT = nT-1:-1:1
%                                 if iT>1
                                    z(:,iT) = z(:,iT+1) + dt(iT) * Dxf(:,:,iT)'*z(:,iT+1) + forcing;
%                                 end
                                JacRaiseF(iX,:) = JacRaiseF(iX,:) + (dt(iT)*Dwf(:,:,iT)'*z(:,iT+1))';
                                if N_alpha > 0
                                    JacRaiseA(iRes,idxAL) = JacRaiseA(iRes,idxAL) + (dt(iT)*Daf(:,:,iT)'*z(:,iT+1))';
                                end 
                            end  
                            if ~problem.fixed_x0  
                                JacRatioI(iX,idxIC) = z(:,1)';  
                            end
                        end
                    end
                end
            end
            if step_finalization && plotPeriod_fun(numIter)
                qq{iS}.y = y;
                qq{iS}.x = x;
            end
%             end
            if q.eval_diff
                iResStart = iResStart + nTeval*nY;
            end
        end
    end    
    
    for ipen = 1:length(pen_handlers)
        pen_handlers{ipen}.out = pen_handlers{ipen}.evaluate(paramsF,paramsG,Alpha,IC,compRes,compGrad,compJac,pen_handlers{ipen},numIter,modelclass);
    end
    
    %ERRdiff = .5*normFactDiff*ERRdiff;
    ERRdiff = normFactDiff*ERRdiff;    
    ERRdiff_ex = normFactDiff*ERRdiff_ex;    
    ERRdiff_test = normFactDiff_test*ERRdiff_test;        
    if alphaPen > 0
        ERRend = normFactEnd*ERRend;
    end     
    if ratioPen > 0
        ERRratio = normFactRatio*ERRratio;
    end
    if alphaRaise > 0
        ERRraise = normFactRaise*sum((XRaise/nS/nT-targetRaise).^2);
        if compRes || compGradRaise           
            ResRaise = XRaise/nS/nT-targetRaise;
        end
    end   
     
    if isempty(idxVarRatio)
        ERRratio = 0;
    end
        
    E = diffPen^2     * ERRdiff     ...
      + ratioPen^2    * ERRratio    ...
      + alphaPen^2    * ERRend      ...
      + alphaRaise^2  * ERRraise    ;
    
    for ipen = 1:length(pen_handlers)
        for ipenloc=1:pen_handlers{ipen}.num_pen
            if pen_handlers{ipen}.coef(ipenloc) > 0
                E = E + pen_handlers{ipen}.coef(ipenloc)^2 * pen_handlers{ipen}.norm_pen(ipenloc) * pen_handlers{ipen}.out.errors(ipenloc);
            end
        end
    end
    
    iOut = 1;    
    if compErr
        varargout{iOut} = E;
        iOut = iOut+1;
    end
    if compGrad  
        if alphaRaise > 0 && compGradRaise
            DF = DF + 2*alphaRaise^2*normFactRaise*JacRaiseF'*ResRaise;
        end
        for ipen = 1:length(pen_handlers)
            for ipenloc = 1:pen_handlers{ipen}.num_pen
                if pen_handlers{ipen}.coef(ipenloc) > 0
                    if pen_handlers{ipen}.dependence_f(ipenloc)
                        DF = DF + pen_handlers{ipen}.coef(ipenloc)^2 * pen_handlers{ipen}.norm_pen(ipenloc) * pen_handlers{ipen}.out.gradF{ipenloc};
                    end
                    if pen_handlers{ipen}.dependence_g(ipenloc)
                        DG = DG + pen_handlers{ipen}.coef(ipenloc)^2 * pen_handlers{ipen}.norm_pen(ipenloc) * pen_handlers{ipen}.out.gradG{ipenloc};
                    end
                    if pen_handlers{ipen}.dependence_a(ipenloc)
                        DA = DA + pen_handlers{ipen}.coef(ipenloc)^2 * pen_handlers{ipen}.norm_pen(ipenloc) * pen_handlers{ipen}.out.gradA{ipenloc};
                    end
                    if pen_handlers{ipen}.dependence_i(ipenloc)
                        DI = DI + pen_handlers{ipen}.coef(ipenloc)^2 * pen_handlers{ipen}.norm_pen(ipenloc) * pen_handlers{ipen}.out.gradI{ipenloc};
                    end
                end
            end
        end
        
        if ~checkGrad
            varargout{iOut} = [DF;DG;DI;DA];
            iOut = iOut+1;
        end
    end
    if compRes
        Residuals = [];
        if diffPen > 0
            Residuals = [Residuals;sqrt(2*normFactDiff)*diffPen*Res];
        end
        if alphaPen > 0
            Residuals = [Residuals;sqrt(2*normFactEnd)*alphaPen*ResEnd];
        end
        if alphaRaise > 0
            Residuals = [Residuals;sqrt(2*normFactRaise)*alphaRaise*ResRaise];
        end
        if ratioPen > 0
            Residuals = [Residuals;sqrt(2*normFactRatio)*ratioPen*ResRatio];
        end
        for ipen = 1:length(pen_handlers)
            for ipenloc = 1:pen_handlers{ipen}.num_pen
                if pen_handlers{ipen}.coef(ipenloc) > 0
                    Residuals = [Residuals; sqrt(2*pen_handlers{ipen}.norm_pen(ipenloc)) * pen_handlers{ipen}.coef(ipenloc) * pen_handlers{ipen}.out.res{ipenloc}];
                end
            end
        end
        
        varargout{iOut} = Residuals;
        iOut = iOut+1;
    end
    if compJac
        Jac = [];
        if diffPen > 0
            Jac = [Jac;sqrt(2*normFactDiff)*diffPen*[JacF JacG JacI JacA]];
        end
        if alphaPen > 0
            Jac = [Jac; sqrt(2*normFactEnd)*alphaPen*[JacEndF JacEndG JacEndI JacEndA]];
        end
        if alphaRaise > 0
            Jac = [Jac; sqrt(2*normFactRaise)*alphaRaise*[JacRaiseF JacRaiseG JacRaiseI JacRaiseA]];
        end
        if ratioPen > 0
            Jac = [Jac; sqrt(2*normFactRatio)*ratioPen*[JacRatioF JacRatioG JacRatioI JacRatioA]];
        end
        for ipen = 1:length(pen_handlers)
            for ipenloc = 1:pen_handlers{ipen}.num_pen
                if pen_handlers{ipen}.coef(ipenloc) > 0
                    residx = size(Jac,1);
                    varidx = 0;
                    varsize = nwF;
                    if pen_handlers{ipen}.dependence_f(ipenloc)
                        Jac(residx+1:residx+pen_handlers{ipen}.size_pen(ipenloc),varidx+1:varidx+varsize) = ...
                            sqrt(2*pen_handlers{ipen}.norm_pen(ipenloc)) * pen_handlers{ipen}.coef(ipenloc) * pen_handlers{ipen}.out.jacF{ipenloc};
                    end
                    varidx = varidx + varsize;
                    varsize = nwG;
                    if pen_handlers{ipen}.dependence_g(ipenloc)
                        Jac(residx+1:residx+pen_handlers{ipen}.size_pen(ipenloc),varidx+1:varidx+varsize) = ...
                            sqrt(2*pen_handlers{ipen}.norm_pen(ipenloc)) * pen_handlers{ipen}.coef(ipenloc) * pen_handlers{ipen}.out.jacG{ipenloc};
                    end
                    varidx = varidx + varsize;
                    varsize = nIC;
                    if pen_handlers{ipen}.dependence_i(ipenloc)
                        Jac(residx+1:residx+pen_handlers{ipen}.size_pen(ipenloc),varidx+1:varidx+varsize) = ...
                            sqrt(2*pen_handlers{ipen}.norm_pen(ipenloc)) * pen_handlers{ipen}.coef(ipenloc) * pen_handlers{ipen}.out.jacI{ipenloc};
                    end
                    varidx = varidx + varsize;
                    varsize = nAL;
                    if pen_handlers{ipen}.dependence_a(ipenloc)
                        Jac(residx+1:residx+pen_handlers{ipen}.size_pen(ipenloc),varidx+1:varidx+varsize) = ...
                            sqrt(2*pen_handlers{ipen}.norm_pen(ipenloc)) * pen_handlers{ipen}.coef(ipenloc) * pen_handlers{ipen}.out.jacA{ipenloc};
                    end
                end
            end
        end
        
        varargout{iOut} = Jac;
        iOut = iOut+1;
    end
    
    
    if checkGrad
        DE = [DF;DG;DI;DA];
        DEbis = Jac'*Residuals; 
        [DE DEbis];
        figure();
        plot(DE,'-','linewidth',2);
        hold on
        plot(DEbis,'--','linewidth',2);
        title(num2str(norm(DE-DEbis)/norm(DE)));
        pause()
        close()
    end

          
    if step_finalization
        E_n = [ERRdiff     ; ...
               ERRdiff_ex  ; ...
               ERRratio    ; ...
               ERRend      ; ...
               ERRraise    ];
        Wgt = [diffPen     ; ...
               diffPen     ; ...
               ratioPen    ; ...
               alphaPen    ; ...
               alphaRaise  ];      
        Leg = {'ERRdiff'        , ...
               'ERRdiff (ex)'   , ...
               'ERRratio'       , ...
               'ERRend'         , ...
               'ERRraise'       };
        Act = [    diffPen_act , ...
                   is_y_noised , ...
                  ratioPen_act , ...
                  alphaPen_act , ...
                alphaRaise_act ];
            
        for ipen = 1:length(pen_handlers)
            for ipenloc=1:pen_handlers{ipen}.num_pen
                if pen_handlers{ipen}.coef(ipenloc) > 0
                    E_n = [E_n; pen_handlers{ipen}.norm_pen(ipenloc) * pen_handlers{ipen}.out.errors(ipenloc)];
                else
                    E_n = [E_n; 0];
                end
                Wgt = [Wgt; pen_handlers{ipen}.coef(ipenloc)];
                Leg{length(Leg)+1} = pen_handlers{ipen}.names{ipenloc};
                Act = [Act, pen_handlers{ipen}.pen_active(ipenloc)];
            end
        end
        EE_t = [EE_t E];
        EE_n = [EE_n E_n];
        EE_w = [EE_w E_n.*Wgt.^2];
        EE_n_test = [EE_n_test ERRdiff_test];
        EE_w_test = [EE_w_test diffPen^2*ERRdiff_test];
        
        xmean_mean   = xmean_mean/nS;  
        xmeanL2_mean = xmeanL2_mean/nS;  
        
        if plotPeriod_fun(numIter) && do_plot
            if numIter>0
                format_legend = @(s,v) sprintf('%s (%1.2e)',s,v);
                legend_labels{1} = format_legend('ERR',sqrt(EE_t(end)));
                for iE = 1:length(Act)
                    if Act(iE)
                        legend_labels = [legend_labels {format_legend(Leg{iE},sqrt(EE_w(iE,end)))}];
                    end
                end
                legend_labels = [legend_labels {format_legend('EFFdiff (test)',sqrt(EE_w_test(end)))}];

                subplot(axHistory_1);
                hold off
                loglog(sqrt(EE_t)','k--','linewidth',2);
                hold on
                loglog(sqrt(EE_w(Act==1,:))','linewidth',2);
                ax = gca;
                ax.ColorOrderIndex = 1;
                loglog(sqrt(EE_w_test)','--','linewidth',2);
                % legend([{'ERR'} Leg{Act==1} {'EFFdiff (test)'}] , ...
                       % 'location','southwest');
                legend(legend_labels,'location','southwest');
                sTitle = strcat('iter=', num2str(numIter));
                if saveNetwork
                    sTitle = strcat(sTitle,' - ',strrep(CurrentDirName,'_','\_'));
                end
                title(sTitle);
                grid on

                if ~problem.samples_variability 
                    subplot(axHistory_2);
                    hold off
                    loglog(sqrt(EE_n(Act==1,:))','linewidth',2);
                    hold on
                    ax = gca;
                    ax.ColorOrderIndex = 1;
                    loglog(sqrt(EE_n_test)','--','linewidth',2);
                    grid on
                end
            end

            if problem.samples_variability 
                subplot(axHistory_2);
                
                if isstruct(normalization.alpha_norm)
                    A_a = (normalization.alpha_norm.max-normalization.alpha_norm.min)/2;
                    B_a = (normalization.alpha_norm.max+normalization.alpha_norm.min)/2;
                else
                    A_a = normalization.alpha_norm;
                    B_a = zeros(N_alpha,1);
                end   
    
                if N_alpha == 1 && isfield(tr{1}, 'alpha')
                    if size(tr{1}.alpha,1) == 1
                        if isfield(normalization,'alpha_to_alpha')
                            npt = 1e2;
                            aa = linspace(0,1,npt);
                            for ia = 1:npt
                                aa_tilde(ia) = normalization.alpha_to_alpha(aa(ia));
                            end
                            plot(aa_tilde,aa,'k-')
                            hold on
                        end
                        for iS = 1:nS
                            hhh = plot(Alpha(1,iS).*A_a + B_a, tr{iS}.alpha,'o');
                            set(hhh, 'MarkerFaceColor', get(hhh,'Color')); 
                            hold on                            
                        end
                        xlabel('learned')
                        ylabel('real')
                        al_min_plot = min(Alpha.*A_a + B_a);
                        al_max_plot = max(Alpha.*A_a + B_a);
                        if al_min_plot == al_max_plot
                            al_min_plot = al_min_plot-1;
                            al_max_plot = al_max_plot+1;
                        end
                        axis([al_min_plot al_max_plot -inf inf])
                        hold off  
                    end
                end
                if N_alpha == 0 && N == nY+1 && size(tr{1}.alpha,1) == 1
                    for iS = 1:nS
                        hhh = plot(IC(end,iS), tr{iS}.alpha,'o');
                        set(hhh, 'MarkerFaceColor', get(hhh,'Color')); 
                        hold on                            
                    end
                    xlabel('learned')
                    ylabel('real')
                    hold off  
                end
                if N_alpha == 2
                    for iS = 1:nS
                        hhh = plot(Alpha(1,iS).*A_a + B_a, Alpha(2,iS).*A_a + B_a,'o');
                        set(hhh, 'MarkerFaceColor', get(hhh,'Color')); 
                        hold on                            
                    end
                    xlabel('\alpha_1')
                    ylabel('\alpha_2')
                    hold off  
                end
            end

            if diffPen > 0 || alphaPen > 0 || alphaRaise > 0 || ratioPen > 0
                %figure(FigMain);
                lineWidth = 1;
                
                for jK = 1:nSOut+nSOut_test
                    
                    if jK <= nSOut
                        subplot(ax_train(jK));
                        q = tr{jK};
                    else
                        subplot(ax_test(jK-nSOut));
                        q = ts{jK-nSOut};
                    end
                    idxPlot = 1:q.nT;
                    hold off
                    plot(q.tt(idxPlot),q.yy(:,idxPlot),'--','linewidth',lineWidth);
                    hold on
                    ax = gca;
                    ax.ColorOrderIndex = 1;
                    plot(q.tt(idxPlot),qq{jK}.y(:,idxPlot),'-','linewidth',lineWidth);

                    if nU > 0 && plot_showU
                        if normalizeplot
                            plot(q.tt,-1 + 2*(q.uu(:,:) - uMin)./(uMax - uMin),'-');
                        else
                            plot(q.tt,q.uu(:,:),'-');
                        end
                    end
                                        
%                     if ~useG
%                         plot(q.tt(idxPlot),qq{jK}.x(1:nY,idxPlot),'-','linewidth',lineWidth);
%                     end

%                     plot(tt(idxPlot),xAll(modelclass.x_idx_plot,idxPlot,jK)./abs(xmean_mean(modelclass.x_idx_plot))*xNormPlot,'-','linewidth',lineWidth);
                    if ~isempty(modelclass.x_idx_plot) && modelclass.plot_x
                        plot(q.tt(idxPlot),qq{jK}.x(modelclass.x_idx_plot,idxPlot)*xNormPlot,'-','linewidth',lineWidth);
                    end
                    if N_alpha > 0
                        plot([q.tt(1) q.tt(end)],Alpha(:,jK)*[1 1],':');
                    end
                    axis([0 q.Tmax yBounds]);
%                     if normalizeplot
                        set(gca,'YTickLabel',[]);
%                     end
                    set(gca,'XTickLabel',[]);
                    
                    if jK > nSOut                              
                        set(gca,'Color',[1 1 1]*.9);
                    end
                end
%             
%                 for jK = 1:nSOut
% %                     subplot(nRows,nCols,getIdxPlot(jK));
%                     subplot(ax_train(jK));
%                     hold off
%                     plot(tt(idxPlot),yy(:,idxPlot,jK),'--','linewidth',lineWidth);
%                     hold on
%                     ax = gca;
%                     ax.ColorOrderIndex = 1;
%                     plot(tt(idxPlot),yAll(:,idxPlot,jK),'-','linewidth',lineWidth);
%                     if nU > 0
%                         plot(tt,-1 + 2*(uu(:,:,jK) - uMin)./(uMax - uMin),'-');
%                     end
%                     if ~useG
%                         plot(tt(idxPlot),xAll(1:nY,idxPlot,jK),'-','linewidth',lineWidth);
%                     end
% %                     plot(tt(idxPlot),xAll(modelclass.x_idx_plot,idxPlot,jK)./abs(xmean_mean(modelclass.x_idx_plot))*xNormPlot,'-','linewidth',lineWidth);
%                     plot(tt(idxPlot),xAll(modelclass.x_idx_plot,idxPlot,jK)*xNormPlot,'-','linewidth',lineWidth);
%                     axis([0 Tmax yBounds]);
%                     set(gca,'YTickLabel',[]);
%                     set(gca,'XTickLabel',[]);
%                 end
% 
%                 for jK = 1:nSOut_test
% %                     subplot(nRows,nCols,getIdxPlot(nSOut+jK));
%                     subplot(ax_test(jK));
%                     hold off
%                     plot(tt(idxPlot),yy_test(:,idxPlot,jK),'--','linewidth',lineWidth);
%                     hold on
%                     ax = gca;
%                     ax.ColorOrderIndex = 1;
%                     plot(tt(idxPlot),yAll(:,idxPlot,nS+jK),'-','linewidth',lineWidth);
%                     if nU > 0
%                         plot(tt, -1 + 2*(uu_test(:,:,jK) - uMin)./(uMax - uMin),'-');
%                     end
%                     if ~useG
%                         plot(tt(idxPlot),xAll(1:nY,idxPlot,nS+jK),'-','linewidth',lineWidth);
%                     end
% %                     plot(tt(idxPlot),xAll(modelclass.x_idx_plot,idxPlot,nS+jK)./abs(xmean_mean(modelclass.x_idx_plot))*xNormPlot,'-','linewidth',lineWidth);
%                     plot(tt(idxPlot),xAll(modelclass.x_idx_plot,idxPlot,nS+jK)*xNormPlot,'-','linewidth',lineWidth);
%                     axis([0 Tmax yBounds]);            
%                     set(gca,'Color',[1 1 1]*.9);
%                     set(gca,'YTickLabel',[]);
%                     set(gca,'XTickLabel',[]);
%                 end
            end


%             if (nY == 1 && nU == 0)
%                 subplot(nRows,nCols,idxANNplot)
%                 vVals = -.2:.01:1.5;
%                 dvVals = zeros(1,length(vVals));
%                 for j = 1:length(vVals)
%                     [dvVals(j),~,~] = EvaluateANN(numnF,wF,thetaF,vVals(j),f,BetaOutput);
%                 end
%                 plot(vVals,dvVals)
%             end
            
    %         if nS == 1
    %             subplot(nRows,nCols,4)
    %             vVals = -.2:.01:1.5;
    %             dvVals = zeros(1,length(vVals));
    %             for j = 1:length(vVals)
    %                 [dvVals(j),~,~] = EvaluateANN(numnF,wF,thetaF,vVals(j),f,BetaOutput);
    %             end
    %             plot(vVals,dvVals)
    %         end

    %         figure(FigANN)
    %         subplot(2,1,1)
    %         VisualizeANN(numnF,wF,thetaF)
    %         subplot(2,1,2)
    %         VisualizeANN(numnG,wG,thetaG)

            pause(1e-16);
        end
        
        
        if modelclass.renormalize_allows
            for iX = modelclass.idx_hidden_states
                normFact = 1;
                if xmean_mean(iX) < 0
                    normFact = -1;
                end
                if xmeanL2_mean(iX) ~= 0 && (xmeanL2_mean(iX) < xmean_liminf || xmeanL2_mean(iX) > xmean_limsup)
                   normFact = normFact * xmean_target / xmeanL2_mean(iX); % if negative do both renormalization and sign swap
                end
                if normFact~=1
                   fprintf('Renormalization of %d-th component by a factor %f...\n',iX,normFact);
                   [paramsF,paramsG] = modelclass.renormalize(normFact,iX,paramsF,paramsG);                                 
                    if ~problem.fixed_x0
                       IC(iX,:) = normFact * IC(iX,:);
                    end
                end
            end
        end
                
        idxCurr = 0;
        ParamsOut(idxCurr+1:idxCurr+nwF)      =  paramsF  ; idxCurr = idxCurr+nwF;
        if useG                                       
            ParamsOut(idxCurr+1:idxCurr+nwG)  =  paramsG  ; idxCurr = idxCurr+nwG;
        end                                                          
        if ~problem.fixed_x0     
            IC = reshape(IC,nIC,1);                               
            ParamsOut(idxCurr+1:idxCurr+nIC)  =  IC       ; idxCurr = idxCurr+nIC;
        end
        if N_alpha > 0
            Alpha_vec = reshape(Alpha,nAL,1);
            ParamsOut(idxCurr+1:idxCurr+nAL)  =  Alpha_vec; idxCurr = idxCurr+nAL;
        end
        
        ParamsOut = ParamsOut';
        varargout{1} = ParamsOut;
                
        if snapPeriod_fun(numIter) && saveNetwork && saveHistory
            fileCurr = sprintf('net_%06d.mat',numIter);
            fileCurrFull = strcat(CurrentDir,'/networks/',fileCurr);
            save(fileCurrFull,'EE_t','EE_n','EE_w','EE_n_test','EE_w_test','Wgt','Leg','Act');            
            if isfield(modelclass,'save_to_file_paramsF')
                modelclass.save_to_file_paramsF(fileCurrFull,paramsF);
            end
            save(fileCurrFull,'paramsF','-append');
            if useG
                if isfield(modelclass,'save_to_file_paramsF')
                    modelclass.save_to_file_paramsG(fileCurrFull,paramsG);
                end
                save(fileCurrFull,'paramsG','-append');
            end  
            if ~problem.fixed_x0
                save(fileCurrFull,'IC','-append');
            end
            if N_alpha > 0
                save(fileCurrFull,'Alpha','-append');
            end
        end
        
        if savePeriod_fun(numIter) && saveNetwork
            save(fileNet,'EE_t','EE_n','EE_w','EE_n_test','EE_w_test','Wgt','Leg','Act');
            if isfield(modelclass,'save_to_file_paramsF')
                modelclass.save_to_file_paramsF(fileNet,paramsF);
            end
            save(fileNet,'paramsF','-append');
            if useG
                if isfield(modelclass,'save_to_file_paramsF')
                    modelclass.save_to_file_paramsG(fileNet,paramsG);
                end
                save(fileNet,'paramsG','-append');
            end 
            if ~problem.fixed_x0
                save(fileNet,'IC','-append');
            end
            if N_alpha > 0
                save(fileNet,'Alpha','-append');
            end
        end
        
%         if mod(numIter,backupPeriod) == 0 && saveNetwork
%             copyfile(OutputFile,char(strcat(extractBefore(OutputFile,'.mat'),'_',num2str(numIter),'.mat')));
%         end
        numIter = numIter+1;
        if saveNetwork
            diary off
            diary(DiaryFile)
        end
    end
end

end