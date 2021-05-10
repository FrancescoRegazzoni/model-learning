clear

% problem = problem_get('testcase','testcase_1var_linear.ini');
problem = problem_get('testcase','testcase_1var_exp.ini');
% problem = problem_get('testcase','testcase_1var_expstab.ini');
% problem = problem_get('testcase','testcase_1var_bistab.ini');

% problem = problem_get('testcase','testcase_1var_expstab_osc.ini');
% problem = problem_get('testcase','testcase_1var_expstab_2a.ini');

% problem = problem_get('testcase','testcase_2var_exp.ini');

HFmod = problem.get_model(problem);
model_check_derivatives(HFmod);

db = database_read('db.txt',problem);
db(db.noise == 0,:).noise = 0*db(db.noise == 0,:).noise + 1e-16;
%db.type_comm = db.type+1
db

%model_show_example(HFmod)

% ANN MOR
% model_learn('opt_testcase_1var.ini')
% model_learn('opt_testcase_2var.ini')
% model_learn('opt_testcase_1var_expstab_2a.ini')
%% Dataset loading
samples = dataset_get(struct('problem',problem,'type','file','source','samples.mat;1:20'));

clear opt_train_plot
opt_train_plot.dummy = 0;
% opt_train_plot.automargin = 0;
% opt_train_plot.figure_size = [500 200 1000 600];
opt_train_plot.show_x_ticks = -1;
opt_train_plot.show_y_ticks = -1;
opt_train_plot.y_label = 'x';
opt_train_plot.num_rows = 4;
opt_train_plot.num_cols = 5;
dataset_plot(samples,problem,opt_train_plot)
%%
type = 0; % 0 = free; 1 = alpha constr.;  2 = alpha penal; 3 = a + df/da penal
noise = 1e-16;
idx = 1;

db_filtered = db(db.noise == noise & db.type == type,:)
ANNmod = read_model_fromfile(problem,db_filtered.dir{idx},struct('show_history',1));
% model_check_derivatives(ANNmod);
model_alpha_plot(ANNmod);
% axis([0 1.1 0 1])
metamodels_comparison(HFmod,ANNmod);

% iS = 4;
% figure()
% plot(samples{iS}.tt,samples{iS}.yy,'k')
% hold on
% plot(samples{iS}.tt,samples{iS}.yy + noise./sqrt(ANNmod.dt).*randn(size(samples{iS}.yy)))
% axis([0 1 .5 2])

%%
figure();
opt_solve.do_plot = 1;
test_solve.tt = [0 problem.T];
ANNmod_part = metamodel_particularize(ANNmod);
out = model_solve(test_solve,ANNmod_part,opt_solve);
%%
da_test(ANNmod,1e-3);
%%
da_test(HFmod,1e-3);
%%
clear opt_dap
opt_dap.mod_HF = HFmod;
opt_dap.fraction_obs = .5;
opt_dap.obs_err = 0;
opt_dap.da_obs_err = 1e-3;
opt_dap.pause_each_test = 0;
opt_dap.do_plot = 1;
% obs_err = 10.^-(0:2:16);  
% obs_err = 10.^-(1:6);
% obs_err = 5e-2;
% HFmod.dt = 1e-3;
% ANNmod.dt = 1e-2;
da_estimate_predict(problem,ANNmod,opt_dap);
% da_estimate_predict(ANNmod,ANNmod,obs_err,struct('da_obs_err',obs_err,'pause_each_test',0,'do_plot',0,'fraction_obs',fraction_obs));
% da_estimate_predict(HFmod,ANNmod,obs_err,struct('da_obs_err',obs_err,'da_mod_err',0,'pause_each_test',0,'do_plot',0,'reference_qty','obs_err','fraction_obs',fraction_obs));
% da_estimate_predict(HFmod,ANNmod,obs_err,struct('da_obs_err',max(obs_err,1e-5),'da_mod_err',1e-6,'pause_each_test',0,'do_plot',0,'reference_qty','da_obs_err'));
%%
subplot(1,3,1)
axis([1e-6 1e-1 1e-6 1e0])
subplot(1,3,2)
axis([1e-6 1e-1 1e-6 1e0])
subplot(1,3,3)
axis([1e-6 1e-1 1e-6 1e3])

%% 
fraction_obs = .5;
type = 1;
% obs_err = 1e-10;
obs_err = 10.^-(1:6);

da_estimate_predict(HFmod,db(db.noise > 2e-16 & db.type == type,:),obs_err,struct('da_mod_err',0,'pause_each_test',1,'do_plot',0,'reference_qty','noise_train','fraction_obs',fraction_obs));

%%
HFmodCN = problem.models('CN').get_model(problem);
model_check_derivatives(HFmodCN);
da_test(HFmodCN,1e-2);