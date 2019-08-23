clear
do_save = 1;
%%
problem = problem_get('transmission_line_circuit','NTL1.ini');
problem.goto_dir()
%% 
dataset_def.problem = problem;
dataset_def.type = 'file';

dataset_def.source = 'samples_x_step.mat;1:5|samples_x_rnd.mat;1:20';
ds_train = dataset_get(dataset_def);
%%
clear opt_train_plot
opt_train_plot.automargin = 0;
opt_train_plot.size_opt.left_margin   = 40;
opt_train_plot.size_opt.right_margin  = 40;
opt_train_plot.size_opt.bottom_margin = 40;
opt_train_plot.size_opt.margin_plot_h = 10;
opt_train_plot.size_opt.margin_plot_v = 5;
opt_train_plot.filename = 'fig/NTL1_trainset';
opt_train_plot.linewidth_u = 1.5;
opt_train_plot.linewidth_y = 1.5;
opt_train_plot.show_x_ticks = -1;

opt_train_plot.do_save = do_save;
if do_save
    create_directory_if_not_found('fig')
end
dataset_plot(ds_train,problem,opt_train_plot)