clear
%% options
do_save = 1;
print_labels = 0;
filename = 'fig/pendulum_test';
horizontaltiling = 1;

cmap = get(0, 'DefaultAxesColorOrder');
col_exact       = cmap(1,:);
col_y           = cmap(2,:);
col_u1          = cmap(3,:);
col_u2          = cmap(4,:);
sty_exact       = '-';
sty_y           = '--';
sty_u1          = ':';
sty_u2          = ':';

%% tests definition
iTest = 1;
tests{iTest}.u1 = @(t) sin(1*t).*cos(1.3*t);
tests{iTest}.u2 = @(t) cos(1.8*t).*sin(1*t);
tests{iTest}.ylim = [-1.1 1.1];
iTest = iTest+1;
tests{iTest}.u1 = @(t) floor((20-t)/4)/5;
tests{iTest}.u2 = @(t) floor(t/4)/5;
tests{iTest}.ylim = [-0.1 1.1];
iTest = iTest+1;
tests{iTest}.u1 = @(t) 0*t+.5;
tests{iTest}.u2 = @(t) (mod(t,8)<4);
tests{iTest}.ylim = [-0.1 1.1];
iTest = iTest+1;
tests{iTest}.u1 = @(t) cos(t);
tests{iTest}.u2 = @(t) sin(t);
tests{iTest}.ylim = [-1.1 1.1];
iTest = iTest+1;

Tmax = 20;
    
%% generate tests
nTests = length(tests);
for iTest = 1:nTests
    tests{iTest}.test.tt = [0 Tmax];
    tests{iTest}.test.uu = @(t) [tests{iTest}.u1(t); tests{iTest}.u2(t)];
end

% load problem and HF model
problem = problem_get('pendulum','pendulum.ini');
HFmod = problem.get_model(problem);

opt_solve.save_x = 1;
for iTest = 1:nTests
    tests{iTest}.sol = model_solve(tests{iTest}.test,HFmod,opt_solve);
end

size_opt.plot_legend = 1;
size_opt.left_margin = 50;
size_opt.top_margin = 30;
size_opt.subplot_width = 150;
size_opt.subplot_heigth = 100;
size_opt.margin_plot_h = 30;
if horizontaltiling
    win_sizes = plotting_get_axes(2,nTests,size_opt);
else
    win_sizes = plotting_get_axes(nTests,2,size_opt);
end
figure('units','pixel','position',[100 100 win_sizes.width win_sizes.heigth]);

load_ANNs;
for N = 1:2        
    for iTest = 1:nTests
        outANN = model_solve(tests{iTest}.test,ANNmodels_best{N},opt_solve);
        
        if horizontaltiling
            iplot = (N-1)*nTests+iTest;
        else
            iplot = 2*(iTest-1)+N;
        end
        subplot('Position',win_sizes.axs{iplot}.coord_norm);
        
        plot(tests{iTest}.sol.tt,tests{iTest}.sol.uu(1,:),sty_u1,'linewidth',1,'Color',col_u1); hold on
        plot(tests{iTest}.sol.tt,tests{iTest}.sol.uu(2,:),sty_u2,'linewidth',1,'Color',col_u2); hold on
        plot(tests{iTest}.sol.tt,tests{iTest}.sol.yy,sty_exact,'linewidth',1,'Color',col_exact); hold on
        plot(outANN.tt,outANN.xx(1,:),sty_y,'linewidth',1,'Color',col_y); hold on 
        
        if N == 1
            title(sprintf('Test #%d',iTest),'FontSize',11)
        end
        if iTest == 1
            ylabel(sprintf('n = %d',N),'FontSize',11,'FontWeight','bold')
        end
        axis([0 Tmax  tests{iTest}.ylim])
    end
end

myleg = legend('u_1','u_2','y (HF)','y (ANN)','Orientation','horizontal');
set(myleg,'position', win_sizes.get_legend_coord_norm(myleg.Position), 'units', 'normalized'); 

if (do_save) 
    create_directory_if_not_found('fig')
    print(filename,'-depsc'); 
end