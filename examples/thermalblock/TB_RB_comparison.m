clear

sourceterm_pb = 0;
if sourceterm_pb
    problem = problem_get('thermalblock','TB9_sourceterm.ini');
else
    problem = problem_get('thermalblock','TB9.ini');
end
problem.goto_dir()
HFmod = problem.get_model(problem);
%% loading ANNs
load_ANNs;
%% comparison ANN vs RB
clear optRBeval
optRBeval.show_epsilon = 1;
optRBeval.normalize_time = 1;
optRBeval.compute_time_HF = 1;
optRBeval.dataset_time_HF = 'samples_rnd.mat;1:5';
optRBeval.comparison_models{1}.models = ANNs_ext;
optRBeval.comparison_models{2}.models = ANNs_int;
optRBeval.comparison_model_name{1} = 'ANN (ext)';
optRBeval.comparison_model_name{2} = 'ANN (int)';
optRBeval.N_values = 1:15;
POD_dataset = 'samples_x_const_lhs80.mat;:|samples_x_rnd_A.mat;1:200|samples_x_rnd_B.mat;1:200';
err_dataset = 'samples_rnd.mat;1:50|samples_x_const_lhs20.mat;:|samples_T10.mat;1:10';
ret = RB_evaluation(HFmod,POD_dataset,err_dataset,optRBeval);
%% saving results
% uncomment the following line to overwrite results
% save('RB_comparison_data.mat','optRBeval','ret');