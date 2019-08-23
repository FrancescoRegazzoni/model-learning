clear

do_save = 1;
inset{1} = [.64 .7 -1.5 -1];
inset{2} = [.73 .82 .6 1];
inset{3} = [.87 .93 .3 .8];
inset{4} = [.61 .71 1.2 1.8];

%%
problem = problem_get('wave_oned','wave1D.ini');
problem.goto_dir()
load_ANNs;
%%
dataset_def.problem = problem;
dataset_def.type = 'file';
dataset_def.source = 'samples_rnd.mat;[33,10,44,47]';
ds = dataset_get(dataset_def);
% ds_out = model_solve(ds,ANNmod);
% dataset_plot(ds_out,problem)

clear size_opt
size_opt.subplot_width = 350;
size_opt.subplot_heigth = 160;
size_opt.margin_plot_v = 40;
size_opt.bottom_margin = 40;
size_opt.plot_legend = 1;
win_sizes = plotting_get_axes(2,2,size_opt);
fig_main = figure('units','pixel','position',[100 100 win_sizes.width win_sizes.heigth]);
win_sizes_mono = plotting_get_axes(1,1,struct('subplot_width',120,'subplot_heigth',80));
for iANN = 1:length(ds)
    fig_mono{iANN} = figure('units','pixel','position',[100 100 win_sizes_mono.width win_sizes_mono.heigth]);
    subplot('Position',win_sizes_mono.axs{1}.coord_norm)
end
cmap = get(0, 'DefaultAxesColorOrder');
lwidU = 1;
lwidYHF = 2;
lwidYANN = 1;
nrows = 2;
ncols = 2;
leg = {'u_1','u_2','y (HF)'};
for iANN = 1:length(ANNs)
    leg = [leg {sprintf('y (ANN, n = %d)',ANNs{iANN}.nX)}];
    ds_out = model_solve(ds,ANNs{iANN});
    for k=1:2        
        for i = 1:length(ds)
            if k == 1
                figure(fig_main);
                subplot('Position',win_sizes.axs{i}.coord_norm)
            else
                figure(fig_mono{i});
            end
            if iANN == 1
                plot(ds_out{i}.tt,ds_out{i}.uu,'--','linewidth',lwidU); hold on;
                plHF = plot(ds_out{i}.tt,ds_out{i}.yy_ex,'k-','linewidth',lwidYHF); hold on;
            end
            plot(ds_out{i}.tt,ds_out{i}.yy,'-','linewidth',lwidYANN); hold on;  
            if k== 1
                xlabel('t'); 
            else
                axis(inset{i})
                set(gca,'XTickLabel',[]);
                set(gca,'YTickLabel',[]);
            end
        end
    end
end

for i = 1:length(ds)
    figure(fig_main);
    subplot('Position',win_sizes.axs{i}.coord_norm)
    in = inset{i};
    col = [1 1 1] * .3;
    plot([in(1) in(2)], [in(3) in(3)], '-', 'Color', col)
    plot([in(2) in(2)], [in(3) in(4)], '-', 'Color', col)
    plot([in(2) in(1)], [in(4) in(4)], '-', 'Color', col)
    plot([in(1) in(1)], [in(4) in(3)], '-', 'Color', col)
end

figure(fig_main);
subplot('Position',win_sizes.axs{1}.coord_norm)
myleg = legend(leg,'Orientation','horizontal');    
set(myleg,'position', win_sizes.get_legend_coord_norm(myleg.Position), 'units', 'normalized');    

if do_save
    create_directory_if_not_found('fig')
    figure(fig_main);
    print('fig/wave1D_samples_test_main','-depsc'); 
    for i = 1:length(ds)
        figure(fig_mono{i});
        print(sprintf('fig/wave1D_samples_test_inset_%d',i),'-depsc'); 
    end
end