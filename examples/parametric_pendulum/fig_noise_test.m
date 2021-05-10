clear
do_table = 1;
n_max_tests = 200;

%% error plot (HF)
data = load('fig/noise_test_HF.mat');
if do_table
    do_table_convergence(data, 'fig/pendulum_noise_test_HF', n_max_tests)
end

%% error plot (online noise = offline noise)
data = load('fig/noise_test_offline_online_noise.mat');
if do_table
    do_table_convergence(data, 'fig/pendulum_noise_test_offline_online_noise', n_max_tests)
end

%% functions

function write_table_points(x,y,xlab,ylab,filename)
    writetable(array2table([reshape(ones(size(y,1),size(y,2)).*x,[],1), reshape(y,[],1)],...
                'VariableNames',{xlab,ylab}),...
                filename,'Delimiter','tab')
end

function do_table_convergence(data, filename_base, n_max_tests)
    arit_mean = mean(data.ret.err_pred_rel_list,1);
    geom_mean = exp(mean(log(data.ret.err_pred_rel_list),1));
    writetable(array2table([data.ret.obs_err', arit_mean', geom_mean'],...
        'VariableNames',{'obs_err','arit_mean','geom_mean'}),...
        [filename_base '_errpred.dat'],'Delimiter','tab')
%         sprintf('fig/testcase_1var_noise_%d_errpred.dat',i_noise),'Delimiter','tab')

    write_table_points(data.ret.obs_err,data.ret.err_pred_rel_list(1:n_max_tests,:), ...
        'obs_err','err_pred', ...
        [filename_base '_errpred_points.dat'])
%         sprintf('fig/testcase_1var_noise_%d_errpred_points.dat',i_noise))
end