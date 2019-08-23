clear
do_save = 0;

%% loading problem
problem = problem_get('thermalblock','TB9.ini');
problem.goto_dir()

%% loading data set
dataset_def.problem = problem;
dataset_def.type = 'file';
dataset_def.source = 'samples_x_const_lhs10.mat;1:10|samples_x_rnd_A.mat;1:50';
ds_train = dataset_get(dataset_def);

%% train set plot
clear opt_train_plot
opt_train_plot.automargin = 0;
% opt_train_plot.figure_size = [500 200 1000 600];
opt_train_plot.show_x_ticks = 1;
opt_train_plot.show_y_ticks = -1;
        
opt_train_plot.size_opt.subplot_width = 200;
opt_train_plot.size_opt.subplot_heigth = 140;
opt_train_plot.size_opt.left_margin   = 40;
opt_train_plot.size_opt.right_margin  = 40;
opt_train_plot.size_opt.bottom_margin = 40;
% opt_train_plot.size_opt.margin_plot_h = 10;
% opt_train_plot.size_opt.margin_plot_v = 5;
opt_train_plot.linewidth_u = 1.5;
opt_train_plot.linewidth_y = 1.5;
opt_train_plot.filename = 'fig/TB_trainset';
opt_train_plot.do_save = do_save;

% opt_train_plot.left_margin   = 5e-2;
% opt_train_plot.right_margin  = 5e-2;
% opt_train_plot.bottom_margin = 5e-2;
% opt_train_plot.top_margin    = 5e-3;
% opt_train_plot.margin_plot_h = 1e-2;
% opt_train_plot.margin_plot_v = 3e-2;
opt_train_plot.idx_samples = [1 4 5 10 11 31 13 14 15 16 53 58];
opt_train_plot.num_rows = 3;
opt_train_plot.num_cols = 4;

if do_save
    create_directory_if_not_found('fig')
end

dataset_plot(ds_train,problem,opt_train_plot)

% print('TB_trainset','-depsc'); 

% clear opt_train_plot
% opt_train_plot.automargin = 0;
% opt_train_plot.figure_size = [500 500 800 500];
% opt_train_plot.left_margin   = 1e-1;
% opt_train_plot.right_margin  = 1e-1;
% opt_train_plot.bottom_margin = 1e-1;
% opt_train_plot.top_margin    = 5e-3;
% opt_train_plot.margin_plot_h = 5e-3;
% opt_train_plot.margin_plot_v = 5e-3;
% opt_train_plot.idx_samples = [1 4 5 10 11 31 13 14 15 16 53 58];
% % opt_train_plot.idx_samples = [1 4 5 10 ];
% opt_train_plot.num_rows = 3;
% opt_train_plot.num_cols = 4;
% dataset_plot(train,problem,opt_train_plot)