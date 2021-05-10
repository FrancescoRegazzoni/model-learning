clear
% problem = problem_get('parametric_pendulum','problems/param_pend_th0.ini');
problem = problem_get('parametric_pendulum','problems/param_pend_th0_c.ini');

HFmod = problem.get_model(problem);

%% Tests generation
rng('default')

optGen.do_plot = 1;
optGen.do_save = 1;
optGen.T = 10;
optGen.outFile = 'samples.mat';
dataset_generate_random(HFmod,500,optGen);