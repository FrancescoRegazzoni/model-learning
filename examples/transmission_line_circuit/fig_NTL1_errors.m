clear
do_save = 0;
%% load problem and HF model
problem = problem_get('transmission_line_circuit','NTL1.ini');
problem.goto_dir()
HFmod = problem.get_model(problem);
%% load reduced models
load_ANNs;
%% 
dataset_def.problem = problem;
dataset_def.type = 'file';
% dataset_def.source = 'samples_x_step.mat;41:50|samples_x_rnd.mat;61:100|samples_T10.mat;:|samples_T100.mat;:';
% dataset_def.source = 'samples_x_step.mat;41:50|samples_x_rnd.mat;61:100|samples_T10.mat;:';
% dataset_def.source = 'samples_x_step.mat;41:50|samples_x_rnd.mat;61:100';
% dataset_def.source = 'samples_x_step.mat;1:5|samples_x_rnd.mat;1:20';
% dataset_def.source = 'samples_x_const.mat;1:50|samples_x_rnd.mat;21:100';
dataset_def.source = 'samples_x_step.mat;26:50|samples_x_rnd.mat;21:100';

dataset_comp = dataset_get(dataset_def);

comp{1}.models = ANNs_ext;
comp{2}.models = ANNs_int;

%%
% optTMP.do_plot = 0;
% optTMP.pause_eachtest = 0;
% model_compute_error(comp{2}.models{1}, dataset_comp,optTMP);
%%
outErr_HF = model_compute_error(HFmod, dataset_comp);

for iComp = 1:length(comp)
    for iM = 1:length(comp{iComp}.models)      
        outErr_comparison = model_compute_error(comp{iComp}.models{iM}, dataset_comp);
        comp{iComp}.comp_nX(iM) = comp{iComp}.models{iM}.nX;
        comp{iComp}.comp_err(iM) = outErr_comparison.err_dataset_L2_norm;
        comp{iComp}.comp_tim(iM) = outErr_comparison.time_mean / outErr_HF.time_mean;
    end
end
%% plot: error comparison
figure('units','pixel','outerposition',[200 200 300 400]);

semilogy(comp{1}.comp_nX,comp{1}.comp_err,'o-','linewidth',1.5)
hold on
semilogy(comp{2}.comp_nX,comp{2}.comp_err,'o-','linewidth',1.5)
grid on
axis([1 3 1e-3 1e-1])
legend('output-outside-the-state','output-inside-the-state','location','south')

xlabel('n')
ylabel('error')

if do_save
    create_directory_if_not_found('fig')
    print('fig/NTL1_errors','-depsc'); 
end

%% plot: time comparison 
figure('units','pixel','outerposition',[200 200 300 400]);

semilogy(comp{1}.comp_nX,comp{1}.comp_tim,'o-','linewidth',1.5)
hold on
semilogy(comp{2}.comp_nX,comp{2}.comp_tim,'o-','linewidth',1.5)

grid on
% axis([1 3 1e-3 1e-1])
legend('output-outside-the-state','output-inside-the-state')
xlabel('n')
ylabel('time_{ANN} / time_{HF}')