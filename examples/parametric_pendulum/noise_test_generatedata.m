clear
do_export = 1;

noise_levels = [1e-1,1e-2,1e-3,1e-4,1e-5];
fraction_obs = .4;
obs_err = noise_levels;
n_tests = 1000;

rng(1)

%%
problem = problem_get('parametric_pendulum','problems/param_pend_th0_c.ini');
HFmod = problem.get_model(problem);
db = database_read('db.txt',problem);
%%
% outputs = [];
% for i_noise = 1:length(noise_levels)
%     noise = noise_levels(i_noise);
%     db_filtered = db(db.noise == noise,:);
%     ANNmod = read_model_fromfile(problem,db_filtered.dir{1},struct('show_history',0));  
%         
%     clear da_opt
%     da_opt.mod_HF = HFmod;
%     da_opt.fraction_obs = fraction_obs;
%     da_opt.obs_err = obs_err;
%     da_opt.da_mod_err = 0;
%     da_opt.pause_each_test = 0;
%     da_opt.do_plot = 0;
%     da_opt.n_tests = n_tests;
%     ret = da_estimate_predict(problem,ANNmod,da_opt);
%     if do_export
%         save(sprintf('fig/noise_test_%1.0e.mat',noise),'noise','ret');
%     end
% end
%%
clear da_opt
da_opt.mod_HF = HFmod;
da_opt.fraction_obs = fraction_obs;
da_opt.obs_err = obs_err;
da_opt.da_mod_err = 0;
da_opt.pause_each_test = 0;
da_opt.do_plot = 0;
da_opt.do_plot_error = 0;
da_opt.n_tests = n_tests;

ret = da_estimate_predict(problem,HFmod,da_opt);
if do_export
    save('fig/noise_test_HF.mat','ret');
end
%%
clear da_opt
da_opt.mod_HF = HFmod;
da_opt.fraction_obs = fraction_obs;
da_opt.obs_err = obs_err;
da_opt.da_mod_err = 0;
da_opt.pause_each_test = 0;
da_opt.do_plot = 0;
da_opt.do_plot_error = 0;
da_opt.reference_qty = 'noise_train';
da_opt.n_tests = n_tests;

ret = da_estimate_predict(problem,db,da_opt);
if do_export
    save('fig/noise_test_offline_online_noise.mat','ret');
end
