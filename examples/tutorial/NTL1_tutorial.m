clear

%% Problem and model definition
problem = problem_get('tutorial','NTL1.ini');
HFmod = problem.get_model(problem);

%% Test the model
test_solve.tt = [0 10];
test_solve.uu = @(t) .5*(1+cos(2*pi*t/10));
figure();
model_solve(test_solve,HFmod,struct('do_plot',1));

%% Generate tests
rng('default')

opt_gen.do_plot = 1;
opt_gen.do_save = 1;
opt_gen.optRandomU.time_scale = .02; 

opt_gen.outFile = 'samples_rnd.mat';
dataset_generate_random(HFmod,100,opt_gen);

opt_gen.constant = 1;
opt_gen.wait_init = 1;
opt_gen.wait_init_time_wait = .2;
opt_gen.wait_init_time_raise = 0;
opt_gen.outFile = 'samples_step.mat';
dataset_generate_random(HFmod,50,opt_gen);

%% Dataset loading
dataset_def.problem = problem;
dataset_def.type = 'file';
dataset_def.source = 'samples_step.mat;1:3|samples_rnd.mat;1:8';
train_dataset = dataset_get(dataset_def);
dataset_plot(train_dataset,problem)

%% Model learning
model_learn('NTL1_opt.ini')

%% ANN model loading
learned_model_name = '...'; %TODO: write here the model name (see tutorial)
ANNmod = read_model_fromfile(problem,learned_model_name);
ANNmod.visualize()

%% ANN model test
figure();
output = model_solve(test_solve,ANNmod,struct('do_plot',1));

%% Error evaluation
dataset_def.problem = problem;
dataset_def.type = 'file';
dataset_def.source = 'samples_step.mat;11:50|samples_rnd.mat;21:100';
test_dataset = dataset_get(dataset_def);
model_compute_error(HFmod, test_dataset);
model_compute_error(ANNmod, test_dataset);
