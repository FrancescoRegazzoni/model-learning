clear

do_save = 1;

nrows = 1;
ncols = 2;

%%
problem = problem_get('wave_oned','wave1D.ini');
problem.goto_dir()
HFmod = problem.get_model(problem);

dataset_def.problem = problem;
dataset_def.type = 'file';
dataset_def.source = 'samples_rnd.mat;10';
ds = dataset_get(dataset_def);
test = ds{1};
out = model_solve(test,HFmod,struct('save_x',1,'do_plot',1));
%
figure('units','pixel','outerposition',[500 500 800 300]);
n = size(out.xx,1)/2;
for iT = round(linspace(1,length(out.tt),13))
    subplot(nrows,ncols,1)
    plot((0:n)/n,out.xx(1:n+1,iT),'-','linewidth',2)
    hold on
    subplot(nrows,ncols,2)
    plot((1:n-1)/n,out.xx(n+2:end,iT),'-','linewidth',2)
    hold on
end
subplot(nrows,ncols,1)
xlabel('$x$', 'interpreter', 'Latex')
ylabel('$\psi(x)$', 'interpreter', 'Latex')
subplot(nrows,ncols,2)
xlabel('$x$', 'interpreter', 'Latex')
ylabel('$\frac{d\psi}{dt}(x)$', 'interpreter', 'Latex')
if do_save
    create_directory_if_not_found('fig')
    print('fig/wave1D_solutions','-depsc');
end