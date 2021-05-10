clear
do_save = 0;
do_plot = 1;
do_table = 1;

do_alpha_plot = 1;
do_error_plot = 1;
do_error_plot_HF = 1;
do_error_plot_offline_online = 1;
do_error_plot_original_learned = 1;

noise_levels = [1e-2,1e-3,1e-4,1e-5];
n_tests = 1e3;
figure_size_alpha = [100 100 1300 220];
figure_size_error = [100 100 1300 250];
figure_size_error_HF = [100 100 300 290];
figure_size_error_original_learned = [100 100 550 300];

dt_train = 1e-2;

%% alpha plots
if do_alpha_plot
    problem = problem_get('testcase','testcase_1var_exp.ini');
    db = database_read('db_figures.txt',problem);
    type = 1; %alpha constr.

    if do_plot
        fig_alpha= figure('units','pixel','position',figure_size_error);
    end
    for i_noise = 1:length(noise_levels)
        noise = noise_levels(i_noise);
        db_filtered = db(db.noise == noise & db.type == type,:);
        ANNmod = read_model_fromfile(problem,db_filtered.dir{1},struct('show_history',0));  
        if do_table
            writetable(array2table([ANNmod.alpha_learned', ANNmod.alpha_original'],...
                'VariableNames',{'alpha_learned','alpha_original'}),...
                sprintf('fig/testcase_1var_noise_%d_alphaplot.dat',i_noise),'Delimiter','tab')
        end
        if do_plot
            figure(fig_alpha)
            subplot(1,length(noise_levels),i_noise)
    %         fig_alpha_plot(ANNmod, sprintf('strong_noise_%1.0e.mat',noise), figure_size_alpha, do_save)
            fig_alpha_plot(ANNmod)
    %         title(sprintf('offline noise = %1.0e',noise/sqrt(dt_train)))
    %         title(sprintf('offline noise = 10^{%d} s^{1/2}',log10(noise)))
    %         title(sprintf('$\\sigma_{\\text{offline}} = 10^{%d} \\, \\text{s}^{1/2}$',log10(noise)), 'Interpreter', 'tex')
            title(sprintf('\\sigma_{{offline}} = 10^{%d} {s}^{1/2}',log10(noise)), 'Interpreter', 'tex')
    %         title('{\sigma}_{offline} = 10^{-4} s^{1/2}', 'Interpreter', 'tex')
            ax = gca;
            ax.FontSize = 11;
        end
    end
    if do_save
        print(fig_alpha, 'fig/testcase_1var_exp_alphaplot_noise','-depsc','-painters');
    end
    pause(1e-16)
end

%% error plots
if do_error_plot
    outputs = [];
    if do_plot
        fig_err = figure('units','pixel','position',figure_size_error);
    end
    for i_noise = 1:length(noise_levels)
        noise = noise_levels(i_noise);
        data = load(sprintf('fig/testcase_1var_noise_%1.0e.mat',noise));
        if do_table
            do_table_convergence(data, sprintf('fig/testcase_1var_noise_%d',i_noise))
%             arit_mean = mean(data.ret.err_pred_rel_list,1);
%             geom_mean = exp(mean(log(data.ret.err_pred_rel_list),1));
%             writetable(array2table([data.ret.obs_err', arit_mean', geom_mean'],...
%                 'VariableNames',{'obs_err','arit_mean','geom_mean'}),...
%                 sprintf('fig/testcase_1var_noise_%d_errpred.dat',i_noise),'Delimiter','tab')
%             
%             write_table_points(data.ret.obs_err,data.ret.err_pred_rel_list(1:200,:), ...
%                 'obs_err','err_pred', ...
%                 sprintf('fig/testcase_1var_noise_%d_errpred_points.dat',i_noise))
        end
        if do_plot
    %         fig_err = figure('units','pixel','position',figure_size_error);
            subplot(1,length(noise_levels),i_noise)
            da_estimate_predict_makeplot(data.ret.obs_err, data.ret.err_pred_rel_list(1:n_tests,:))
            axis equal
            axis([1e-6 1 1e-6 1])
    %         title(sprintf('offline noise = %1.0e',noise/sqrt(dt_train)))
    %         xlabel('online noise [s^{1/2}]')
            xlabel('\sigma_{online}  [s^{1/2}]')
            ylabel('prediction error [-]')
            title(sprintf('\\sigma_{{offline}} = 10^{%d} {s}^{1/2}',log10(noise)), 'Interpreter', 'tex')
    %         if do_save
    %             print(modelfig, sprintf('fig/testcase_1var_exp_predictionerror_noise_%1.0e', noise),'-depsc','-painters');
    %         end
            ax = gca;
            ax.FontSize = 11;
        end
    end
    if do_plot
        set(fig_err,'PaperPositionMode','auto')
    end
    pause(1e-16);
    if do_save
        print(fig_err, 'fig/testcase_1var_exp_predictionerror_noise','-depsc','-painters');
    end
    pause(1e-16)
end
%% error plot (HF)
if do_error_plot_HF
    data = load('fig/testcase_1var_HF.mat');
    if do_table
        do_table_convergence(data, 'fig/testcase_1var_HF')
    end
    if do_plot
        fig_err = figure('units','pixel','position',figure_size_error_HF);

        da_estimate_predict_makeplot(data.ret.obs_err, data.ret.err_pred_rel_list(1:n_tests,:))
        title('Original model')
    %     xlabel('online noise [s^{1/2}]')
        xlabel('\sigma_{online}  [s^{1/2}]')
        ylabel('prediction error [-]')
        set(fig_err,'PaperPositionMode','auto')
        axis([1e-5 1e-1 1e-5  1])
        pause(1e-16);
        xticks(10.^(-5:-1))
        yticks(10.^(-5:0))

        if do_save
            print(fig_err, 'fig/testcase_1var_exp_predictionerror_HF','-depsc','-painters');
        end
    end
end
%% error plot (online noise = offline noise)

if do_error_plot_offline_online
    data = load('fig/testcase_1var_offline_online_noise.mat');
    if do_table
        do_table_convergence(data, 'fig/testcase_1var_offline_online_noise')
    end
    if do_plot
        fig_err = figure('units','pixel','position',figure_size_error_HF);

        da_estimate_predict_makeplot(data.ret.obs_err, data.ret.err_pred_rel_list(1:n_tests,:))
        title('Learned model')
    %     xlabel('online noise = offline noise [s^{1/2}]')
        xlabel('\sigma_{online} = \sigma_{offline}  [s^{1/2}]')
        ylabel('prediction error [-]')
        set(fig_err,'PaperPositionMode','auto')
        axis([1e-5 1e-1 1e-5  1])
        pause(1e-16);
        xticks(10.^(-5:-1))
        yticks(10.^(-5:0))

        if do_save
            print(fig_err, 'fig/testcase_1var_exp_offline_online_noise','-depsc','-painters');
        end
    end
end
%%
if do_error_plot_original_learned
    if do_plot
        data_HF = load('fig/testcase_1var_HF.mat');
        data_ANN = load('fig/testcase_1var_offline_online_noise.mat');
        fig_err = figure('units','pixel','position',figure_size_error_original_learned );

        c_ord = get(gca,'colororder');

        loglog(data_ANN.ret.obs_err',mean(data_ANN.ret.err_pred_rel_list(1:n_tests,:),1)','.--','linewidth',1.2, 'MarkerSize',25,'Color',c_ord(1,:)); hold on
        loglog(data_ANN.ret.obs_err',exp(mean(log(data_ANN.ret.err_pred_rel_list(1:n_tests,:)),1))','.--','linewidth',1.2, 'MarkerSize',25,'Color',c_ord(2,:)); hold on
        loglog(data_HF.ret.obs_err',mean(data_HF.ret.err_pred_rel_list(1:n_tests,:),1)','.-','linewidth',1.2, 'MarkerSize',25,'Color',c_ord(1,:)); hold on
        loglog(data_HF.ret.obs_err',exp(mean(log(data_HF.ret.err_pred_rel_list(1:n_tests,:)),1))','.-','linewidth',1.2, 'MarkerSize',25,'Color',c_ord(2,:)); hold on
        grid on
        axis([1e-5 1e-1 1e-5  1])

        title('Learned vs original')
    %     xlabel('online noise [s^{1/2}]')
        xlabel('\sigma_{online} = \sigma_{offline}  [s^{1/2}]')
        ylabel('prediction error [-]')
        set(fig_err,'PaperPositionMode','auto')
        pause(1e-16);
        legend('Learned mod. (ar. mean)',...
               'Learned mod. (geom. mean)',...
               'Original mod. (ar. mean)',...
               'Original mod. (geom. mean)',...
               'location','eastoutside')
        xticks(10.^(-5:-1))
        yticks(10.^(-5:0))
        if do_save
            print(fig_err, 'fig/testcase_1var_exp_predictionerror_HFvsANN','-depsc','-painters');
        end
    end
end

%%
function da_estimate_predict_makeplot(x_values,err_list)    
    loglog(x_values',err_list','k+'); hold on
%     p_max = loglog(x_values',max(err_list,[],1)','bo-','linewidth',1.2);
%     p_min = loglog(x_values',min(err_list,[],1)','bo-','linewidth',1.2);
    p_arit_mean = loglog(x_values',mean(err_list,1)','.-','linewidth',1.2, 'MarkerSize',25);
    p_geom_mean = loglog(x_values',exp(mean(log(err_list),1))','.-','linewidth',1.2, 'MarkerSize',25);
    grid on
    hold off
    pause(1e-16);
end

function write_table_points(x,y,xlab,ylab,filename)
    writetable(array2table([reshape(ones(size(y,1),size(y,2)).*x,[],1), reshape(y,[],1)],...
                'VariableNames',{xlab,ylab}),...
                filename,'Delimiter','tab')
end

function do_table_convergence(data, filename_base)
    arit_mean = mean(data.ret.err_pred_rel_list,1);
    geom_mean = exp(mean(log(data.ret.err_pred_rel_list),1));
    writetable(array2table([data.ret.obs_err', arit_mean', geom_mean'],...
        'VariableNames',{'obs_err','arit_mean','geom_mean'}),...
        [filename_base '_errpred.dat'],'Delimiter','tab')
%         sprintf('fig/testcase_1var_noise_%d_errpred.dat',i_noise),'Delimiter','tab')

    write_table_points(data.ret.obs_err,data.ret.err_pred_rel_list(1:200,:), ...
        'obs_err','err_pred', ...
        [filename_base '_errpred_points.dat'])
%         sprintf('fig/testcase_1var_noise_%d_errpred_points.dat',i_noise))
end