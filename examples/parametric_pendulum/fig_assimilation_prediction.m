clear
noise_level = 1e-3;
alpha = [-.1;.6];
x0 = [.8;.1];
do_table = 1;

rng(6)
%%
problem = problem_get('parametric_pendulum','problems/param_pend_th0_c.ini');
HFmod = problem.get_model(problem);
db = database_read('db.txt',problem);
db_filtered = db(db.noise == noise_level,:);
mod_learned = read_model_fromfile(problem,db_filtered.dir{1},struct('show_history',0));

names_da = {'t','x1_exact','x2_exact','x1_estimated','x2_estimated','alpha1','alpha2','x1_std','x2_std','alpha1_std','alpha2_std'};
names_pred = {'t','x1_predicted','x2_predicted','x1_exact','x2_exact'};

opt_dap.dataset = {model_solve(struct('tt',[0 10]), metamodel_particularize(HFmod,x0,alpha))};
opt_dap.dataset_add_noise = 1;
opt_dap.mod_HF = HFmod;
opt_dap.fraction_obs = .4;
opt_dap.obs_err = noise_level;
opt_dap.da_obs_err = noise_level;
opt_dap.do_plot = 0;
opt_dap.out_last_data_assimilation_prediction = 1;
opt_dap.n_tests = 1;
opt_dap.out_alpha = 1;
opt_dap.do_plot = 1;
ret = da_estimate_predict(problem,mod_learned,opt_dap);

if do_table
    writetable(array2table([ret.out_da.tt', ret.out_da.yy', ret.out_da.xx', ret.out_da.PP'],...
        'VariableNames',names_da),...
        'fig/pendulum_assimilation_prediction_data_assimilation.dat','Delimiter','tab')

    writetable(array2table([ret.out_predict.tt', ret.out_predict.yy', ret.out_predict.yy_ex'],...
        'VariableNames',names_pred),...
    'fig/pendulum_assimilation_prediction_prediction.dat','Delimiter','tab')
end