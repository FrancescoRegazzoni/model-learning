clear

problem = problem_get('thermalblock','TB9.ini');
% problem = problem_get('thermalblock','TB9_sourceterm.ini');

HFmod = problem.get_model(problem);

clear optGen_base

optGen_base.do_plot = 1;
optGen_base.do_save = 1;
optGen_base.append = 0;
optGen_base.save_x_dt = .1;
optGen_base.pause_eachtest = 0;
optGen_base.optRandomU.time_scale = .1; 
%%
rng('default')
optGen = optGen_base;
optGen.T = 1;
optGen.constant = 0;
optGen.save_x = 1;
% generate two separate datasets to contain the size of the file
optGen.outFile = 'samples_x_rnd_A.mat';
dataset_generate_random(HFmod,200,optGen);
optGen.outFile = 'samples_x_rnd_B.mat';
dataset_generate_random(HFmod,200,optGen);
%%
rng('default')
optGen = optGen_base;
optGen.T = 1;
optGen.constant = 0;
optGen.save_x = 1;
optGen.save_x_dt = 1e-2;
optGen.outFile = 'samples_x_detailed_rnd.mat';
dataset_generate_random(HFmod,1,optGen);
%%
rng(1)
optGen = optGen_base;
optGen.T = 1;
optGen.constant = 0;
optGen.save_x = 1;
optGen.outFile = 'samples_rnd.mat';
dataset_generate_random(HFmod,200,optGen);
%%
for n = [5, 10, 20, 40, 80]
    rng('default')
    optGen = optGen_base;
    optGen.T = .4;
    optGen.constant = 1;
    optGen.lhs = 1;
    optGen.save_x = 1;
    optGen.outFile = sprintf('samples_x_const_lhs%d.mat',n);
    dataset_generate_random(HFmod,n,optGen);
end
%%
rng(2)
optGen = optGen_base;
optGen.T = 10;
optGen.constant = 0;
optGen.save_x = 0;
optGen.outFile = 'samples_T10.mat';
dataset_generate_random(HFmod,10,optGen);
