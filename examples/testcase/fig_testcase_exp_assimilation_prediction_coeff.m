clear
dirANN{1} = 'test_int_N1_hlayF3_dof13_2019-01-28_10-29-25';
dirANN{2} = 'test_int_N2_hlayF4_dof26_2019-02-19_18-36-35';

for i = 1:2
    %%
    problem = problem_get('testcase',sprintf('testcase_%dvar_exp.ini',i));
    HFmod = problem.get_model(problem);
    mod_learned = read_model_fromfile(problem,dirANN{i});
    writetable(array2table([mod_learned.alpha_learned', mod_learned.alpha_original'],...
        'VariableNames',{'alpha_learned','alpha_original'}),...
        sprintf('fig/testcase_%dvar_exp_assimilation_prediction_coeff.dat',i),'Delimiter','tab')
    %%
    opt_dap_base.mod_HF = HFmod;
    opt_dap_base.fraction_obs = .5;
    opt_dap_base.obs_err = 0;
    opt_dap_base.da_obs_err = 1e-3;
    opt_dap_base.do_plot = 0;
    %%
    opt_dap = opt_dap_base;
    opt_dap.out_alpha = 1;
    opt_dap.n_tests = 40;
    rng('default')
    ret = da_estimate_predict(problem,mod_learned,opt_dap);
    writetable(array2table([ret.alpha_estimated', ret.alpha_original'],...
        'VariableNames',{'alpha_estimated','alpha_original'}),...
        sprintf('fig/testcase_%dvar_exp_assimilation_prediction_estimated.dat',i),'Delimiter','tab')
    %%
    if i==1
        rng(1)
        names_da = {'t','x_exact','x_estimated','alpha','x_std','alpha_std'};
        names_pred = {'t','x_predicted','x_exact'};
    else
        rng(4)
        names_da = {'t','x1_exact','x2_exact','x1_estimated','x2_estimated','alpha','x1_std','x2_std','alpha_std'};
        names_pred = {'t','x1_predicted','x2_predicted','x1_exact','x2_exact'};
    end
    opt_dap = opt_dap_base;
    opt_dap.out_last_data_assimilation_prediction = 1;
    opt_dap.n_tests = 1;
    opt_dap.out_alpha = 1;
    opt_dap.do_plot = 1;
    ret = da_estimate_predict(problem,mod_learned,opt_dap);

    writetable(array2table([ret.alpha_estimated', ret.alpha_original'],...
        'VariableNames',{'alpha_estimated','alpha_original'}),...
        sprintf('fig/testcase_%dvar_exp_assimilation_prediction_estimated_single.dat',i),'Delimiter','tab')

    writetable(array2table([ret.out_da.tt', ret.out_da.yy', ret.out_da.xx', ret.out_da.PP'],...
        'VariableNames',names_da),...
        sprintf('fig/testcase_%dvar_exp_assimilation_prediction_data_assimilation.dat',i),'Delimiter','tab')

    writetable(array2table([ret.out_predict.tt', ret.out_predict.yy', ret.out_predict.yy_ex'],...
        'VariableNames',names_pred),...
        sprintf('fig/testcase_%dvar_exp_assimilation_prediction_prediction.dat',i),'Delimiter','tab')
end

