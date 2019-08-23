clear
do_save = 1;
%% problem definition
problem = problem_get('pendulum','pendulum.ini');

%% datasets loading
dataset_def.problem = problem;
dataset_def.type = 'file';
dataset_def.source = 'train.mat;:';
train = dataset_get(dataset_def);
dataset_def.source = 'validation.mat;:';
validation = dataset_get(dataset_def);

%% datasets plot
if do_save
    create_directory_if_not_found('fig')
end

opt_plot.show_x_ticks = 1;
opt_plot.show_y_ticks = -1;
opt_plot.uy_same_axis = 1;
opt_plot.linewidth_u = 1.5;
opt_plot.linewidth_y = 1.5;
opt_plot.do_save = do_save;
% opt_plot.do_saveexternallegend = 1;
opt_plot.plot_legend = 1;
opt_plot.size_opt.subplot_width = 120;
opt_plot.size_opt.subplot_heigth = 80;

opt_plot.automargin = 0;
% opt_plot.left_margin   = 1e-1;
% opt_plot.right_margin  = 1e-1;
% opt_plot.top_margin    = 5e-3;
% opt_plot.margin_plot_h = 1e-2;
% opt_plot.margin_plot_v = 2e-2;
% 
% opt_plot.bottom_margin = 2e-2;
% opt_plot.figure_size = [500 200 800 600];
opt_plot.num_rows = 4;
opt_plot.num_cols = 4;
opt_plot.filename = 'fig/samples_trainset';
dataset_plot(train,problem,opt_plot)

% opt_plot.bottom_margin = 2e-1;
% opt_plot.figure_size = [300 500 900 250];
opt_plot.num_rows = 1;
opt_plot.num_cols = 5;
opt_plot.filename = 'fig/samples_validationset';
dataset_plot(validation,problem,opt_plot)
