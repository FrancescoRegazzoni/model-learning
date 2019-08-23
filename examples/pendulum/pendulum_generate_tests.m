clear

% problem definition
problem = problem_get('pendulum','pendulum.ini');
% Model generation
HFmod = problem.get_model(problem);

%% Tests generation
optGen.do_save = 1;

optGen.outFile = 'train.mat';
dataset_generate(HFmod,get_dataset('train'),optGen);

optGen.outFile = 'validation.mat';
dataset_generate(HFmod,get_dataset('validation'),optGen);

optGen.outFile = 'test.mat';
dataset_generate(HFmod,get_dataset('test'),optGen);