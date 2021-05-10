clear

% problem = problem_get('testcase','testcase_1var_linear.ini');
% problem = problem_get('testcase','testcase_1var_exp.ini');
% problem = problem_get('testcase','testcase_1var_expstab.ini');
% problem = problem_get('testcase','testcase_1var_bistab.ini');

% problem = problem_get('testcase','testcase_1var_expstab_osc.ini');
% problem = problem_get('testcase','testcase_1var_expstab_2a.ini');

problem = problem_get('testcase','testcase_2var_exp.ini');

HFmod = problem.get_model(problem);

%% Tests generation
rng('default')

optGen.do_plot = 1;
optGen.do_save = 1;
optGen.T = 1;
optGen.outFile = 'samples.mat';
dataset_generate_random(HFmod,100,optGen);