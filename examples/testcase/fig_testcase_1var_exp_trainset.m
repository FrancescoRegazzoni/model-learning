clear
do_save = 0;
do_table = 1;
%%
problem = problem_get('testcase','testcase_1var_exp.ini');
samples = dataset_get(struct('problem',problem,'type','file','source','samples.mat;1:20'));
%%
if do_save
    opt_train_plot.size_opt.subplot_width = 100;
    opt_train_plot.size_opt.subplot_heigth = 60;
    opt_train_plot.size_opt.left_margin   = 40;
    opt_train_plot.size_opt.right_margin  = 40;
    opt_train_plot.size_opt.bottom_margin = 40;
    % opt_train_plot.size_opt.margin_plot_h = 10;
    % opt_train_plot.size_opt.margin_plot_v = 5;
    opt_train_plot.show_x_ticks = -1;
    opt_train_plot.show_y_ticks = -1;
    opt_train_plot.y_label = '\theta';
    opt_train_plot.num_rows = 4;
    opt_train_plot.num_cols = 5;
    opt_train_plot.grid_on = 1;
    opt_train_plot.linewidth_u = 1.5;
    opt_train_plot.linewidth_y = 1.5;
    opt_train_plot.filename = 'fig/testcase_1var_exp_trainset';
    opt_train_plot.do_save = do_save;
    dataset_plot(samples,problem,opt_train_plot)
end
%%
if do_table
    dataset_export_table(problem, samples, 'fig/testcase_1var_exp_trainset')
end