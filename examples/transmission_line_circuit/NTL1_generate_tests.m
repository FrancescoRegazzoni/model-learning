clear

problem = problem_get('transmission_line_circuit','NTL1.ini');
HFmod = problem.get_model(problem);

%% Tests generation
rng('default')

optGen.do_plot = 1;
optGen.do_save = 1;
optGen.append = 0;
optGen.save_x_dt = .1;
optGen.save_x = 0;
optGen.pause_eachtest = 0;
optGen.optRandomU.time_scale = .02; 
%
optGen.T = 1;
optGen.constant = 0;
optGen.save_x = 1;
optGen.outFile = 'samples_x_rnd.mat';
dataset_generate_random(HFmod,100,optGen);
%
% optGen.T = 1;
% optGen.constant = 1;
% optGen.save_x = 1;
% optGen.outFile = 'samples_x_const.mat';
% dataset_generate_random(HFmod,50,optGen);
%
optGen.T = 1;
optGen.constant = 1;
optGen.save_x = 1;
optGen.wait_init = 1;
optGen.wait_init_time_wait = .2;
optGen.wait_init_time_raise = 0;
optGen.outFile = 'samples_x_step.mat';
dataset_generate_random(HFmod,50,optGen);
%
% optGen.T = 10;
% optGen.constant = 0;
% optGen.save_x = 0;
% optGen.outFile = 'samples_T10.mat';
% dataset_generate_random(HFmod,10,optGen);
%
% optGen.T = 100;
% optGen.constant = 0;
% optGen.save_x = 0;
% optGen.outFile = 'samples_T100.mat';
% dataset_generate_random(HFmod,2,optGen);
%
% optGen.T = 19;
% optGen.constant = 1;
% optGen.closed_loop = 1;
% optGen.closed_loop_toll_x = 2e-1;
% optGen.closed_loop_toll_y = 1e-2;
% optGen.closed_loop_timerelax = 1;
% optGen.closed_loop_timerest = 9;
% optGen.outFile = 'samples_loop_const.mat';
% dataset_generate_random(HFmod,20,optGen);
% optGen.closed_loop = 0;
% %
% optGen.T = 10;
% optGen.constant = 0;
% optGen.closed_loop = 1;
% optGen.closed_loop_toll_x = 2e-1;
% optGen.closed_loop_toll_y = 1e-2;
% optGen.closed_loop_timerelax = 1;
% optGen.closed_loop_timerest = 4;
% optGen.outFile = 'samples_loop.mat';
% dataset_generate_random(HFmod,50,optGen);
% optGen.closed_loop = 0;
