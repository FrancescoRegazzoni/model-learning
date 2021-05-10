clear
% problem = problem_get('parametric_pendulum','problems/param_pend_th0.ini');
problem = problem_get('parametric_pendulum','problems/param_pend_th0_c.ini');
HFmod = problem.get_model(problem);
model_check_derivatives(HFmod);

db = database_read('db.txt',problem);
db(db.noise == 0,:).noise = 0*db(db.noise == 0,:).noise + 1e-16;
db

% model_show_example(HFmod)
%% Dataset loading
samples = dataset_get(struct('problem',problem,'type','file','source','samples.mat;1:20'));

clear opt_train_plot
opt_train_plot.dummy = 0;
opt_train_plot.show_x_ticks = -1;
opt_train_plot.show_y_ticks = -1;
opt_train_plot.y_label = 'x';
opt_train_plot.num_rows = 4;
opt_train_plot.num_cols = 5;
dataset_plot(samples,problem,opt_train_plot)
%%
% model_learn('opt_param_pend_th0_c.ini')
%% Load an ANN model
noise = 1e-4;
idx = 1;

db_filtered = db(db.noise == noise,:)
ANNmod = read_model_fromfile(problem,db_filtered.dir{idx},struct('show_history',1));
model_check_derivatives(ANNmod);
%% Show parameters (learned versus true)
model_alpha_plot(ANNmod);
%% Show fit on training set
model_get_train_fit(ANNmod);
%%
figure();
opt_solve.do_plot = 1;
test_solve.tt = [0 problem.T];
ANNmod_part = metamodel_particularize(ANNmod);
out = model_solve(test_solve,ANNmod_part,opt_solve);
%% Test data-assimilation with HF model
da_test(HFmod,1e-2);
%% Test data-assimilation with ROM
da_test(ANNmod,1e-2);
%% Test estimate-predict framework (single ANN model)

da_opt.mod_HF = HFmod;
da_opt.obs_err = 10.^-(1:8);
% da_opt.obs_err = 10.^-(4);
da_opt.da_obs_err = da_opt.obs_err;
da_opt.da_mod_err = 0;
da_opt.pause_each_test = 0;
da_opt.do_plot = 0;
da_opt.reference_qty = 'obs_err';
da_opt.fraction_obs = .4;

da_estimate_predict(problem,ANNmod,da_opt);

%% Test estimate-predict framework (many ANN models, wrt observation error)
clear da_opt
da_opt.mod_HF = HFmod;
da_opt.obs_err = 10.^-(1:5);
da_opt.fraction_obs = .4;
da_opt.da_mod_err = 0;
da_opt.pause_each_test = 0;
da_opt.do_plot = 0;
da_opt.reference_qty = 'obs_err';

da_estimate_predict(problem,db(:,:),da_opt);
% da_estimate_predict(problem,HFmod,da_opt);

%% Test estimate-predict framework (many ANN models, wrt d.a. model error)

da_opt.mod_HF = HFmod;
da_opt.obs_err = 10.^-(4);
da_opt.da_obs_err = da_opt.obs_err;
da_opt.da_mod_err = 10.^-(14:-2:2);
da_opt.pause_each_test = 0;
da_opt.do_plot = 0;
da_opt.reference_qty = 'da_mod_err';
da_opt.fraction_obs = .4;

da_estimate_predict(problem,ANNmod,da_opt);