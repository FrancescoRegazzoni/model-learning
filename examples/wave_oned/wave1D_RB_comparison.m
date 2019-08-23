clear

problem = problem_get('wave_oned','wave1D.ini');
problem.goto_dir()
HFmod = problem.get_model(problem);
%%
load_ANNs;
%% ANN vs RB comparison
clear optRBeval
optRBeval.show_epsilon = 1;
optRBeval.normalize_time = 1;
optRBeval.compute_time_HF = 1;
optRBeval.dataset_time_HF = 'samples_rnd.mat;1:5';
optRBeval.comparison_models = ANNs;
optRBeval.comparison_model_name = 'ANN';
% optRBeval.comparison_models_times_nskip = 20;
% optRBeval.comparison_models_times_ncompute = 20;
optRBeval.N_values = [1:10 12:2:20 25:5:50 60:10:100 150:50:500];
% optRBeval.N_values = 1;
POD_dataset = 'samples_x_rnd.mat;1:100';
err_dataset = 'samples_rnd.mat;:';
ret = RB_evaluation(HFmod,POD_dataset,err_dataset,optRBeval);

%% viusualization of modes
POD_dataset = 'samples_x_rnd.mat;1:200';
X = build_snapshots_matrix(problem,POD_dataset);

% optPOD.toll = 1e-5;
% optPOD.N = 5;
% optPOD.fixed_N = 1;
optPOD.get_full_V = 1;
optPOD.do_plot = 1;
optPOD.show_epsilon = 1;
[V,outPOD] = POD_projection(X,optPOD);
%
figure()
for i = 1:size(X,2)
    HFmod.make_plot(V(:,i));
    title(sprintf('Mode nr %d; epsilon = %1.2e',i,POD_getepsilon(outPOD.s,i)))
    pause()
end