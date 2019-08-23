clear

do_save = 1;
nSVD = 9;

problem = problem_get('wave_oned','wave1D.ini');
problem.goto_dir()
HFmod = problem.get_model(problem);

%%%%%% RB comparison
% clear optRBeval
% optRBeval.show_epsilon = 1;
% optRBeval.normalize_time = 1;
% optRBeval.compute_time_HF = 1;
% optRBeval.dataset_time_HF = 'samples_rnd.mat;1:5';
% optRBeval.comparison_models = ANNs;
% optRBeval.comparison_model_name = 'ANN';
% % optRBeval.comparison_models_times_nskip = 20;
% % optRBeval.comparison_models_times_ncompute = 20;
% % optRBeval.N_values = [1:10 12:2:20 25:5:50 60:10:100 150:50:500];
% optRBeval.N_values = 1;
% POD_dataset = 'samples_x_rnd.mat;1:100';
% err_dataset = 'samples_rnd.mat;:';
% ret = RB_evaluation(HFmod,POD_dataset,err_dataset,optRBeval);


dataset_def.problem = problem;
dataset_def.type = 'file';
dataset_def.source = 'samples_rnd.mat;:';
ds_test = dataset_get(dataset_def);

load_ANNs;
nmod = length(ANNs);
for i = 1:nmod
    out= model_compute_error(ANNs{i},ds_test);
    err(i) = out.err_dataset_L2_norm;
    Nvars(i) = ANNs{i}.nX;
end

POD_dataset = 'samples_x_rnd.mat;1:100';
X = build_snapshots_matrix(problem,POD_dataset);

optPOD.get_full_V = 1;
[V,outSVD] = POD_projection(X,optPOD);
%%
for iN = 1:nSVD
    ep(iN) = POD_getepsilon(outSVD.s,iN);
end
%%
figure('units','pixel','outerposition',[200 200 500 400]);
semilogy(1:nSVD, ep,'k--','linewidth',1.2)
hold on
semilogy(Nvars, err,'o-','linewidth',1.5)
axis([1 nSVD .5e-3 1])
grid on
legend('\epsilon_n','Test error (ANN)','Location','southwest')
xlabel('n')

if do_save
    create_directory_if_not_found('fig')
    print('fig/wave1D_testerr','-depsc');
end
