clear

problem = problem_get('tutorial','expmod.ini');
HFmod = problem.get_model(problem);

%% Model test
HFmod_particularized = metamodel_particularize(HFmod,1.1,.5);
figure()
model_solve(struct('tt',[0 1]),HFmod_particularized,struct('do_plot',1))

model_show_example(HFmod);
%% Tests generation
rng('default')

opt_gen.do_plot = 1;
opt_gen.do_save = 1;
opt_gen.T = 1;
opt_gen.outFile = 'samples.mat';
dataset_generate_random(HFmod,100,opt_gen);

%% Model learning
model_learn('expmod_opt.ini')

%% Loading learned model
learned_model_name = '...'; %TODO: write here the model name (see tutorial)
ANNmod = read_model_fromfile(problem,learned_model_name);
model_alpha_plot(ANNmod);
model_show_example(ANNmod);

%% Data assimilation test
da_test(HFmod,1e-3);
da_test(ANNmod,1e-3);

%% Data assimilation + prediction
opt_ep.mod_HF = HFmod;
opt_ep.obs_err = 1e-3;
opt_ep.pause_each_test = 1;
opt_ep.do_plot = 1;
da_estimate_predict(problem,ANNmod,opt_ep);

%% Data assimilation + prediction (convergence)
opt_ep.mod_HF = HFmod;
opt_ep.obs_err = 10.^-(1:8);
opt_ep.pause_each_test = 0;
opt_ep.do_plot = 0;
opt_ep.n_tests = 20;
da_estimate_predict(problem,ANNmod,opt_ep);