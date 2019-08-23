clear

do_save = 0;
do_plot_u = 0;
do_inset = 0;

% sty_u = '--';
% sty_HF = 'k-';
% sty_ROM = '-';
% lwidU = 1;
% lwidYHF = 2;
% lwidYANN = 1;

sty_u = '--';
sty_HF = 'k-';
sty_ROM = '--';
lwidU = 1;
lwidYHF = 1.5;
lwidYANN = 1.5;
nrows = 1;
ncols = 2;
%%
problem = problem_get('thermalblock','TB9.ini');
problem.goto_dir()
HFmod = problem.get_model(problem);

dataset_def.problem = problem;
dataset_def.type = 'file';
dataset_def.source = 'samples_rnd.mat;[1,4]';
ds = dataset_get(dataset_def);
% ds_out = model_solve(ds,ANNmod);
% dataset_plot(ds_out,problem)

%%    
sourceterm_pb = 0;
load_ANNs;
%%
POD_dataset = 'samples_x_const_lhs80.mat;:|samples_x_rnd_A.mat;1:200|samples_x_rnd_B.mat;1:200';
X = build_snapshots_matrix(problem,POD_dataset);
%% get RB models
optPOD.get_full_V = 1;
[Vfull,outputPOD] = POD_projection(X,optPOD);
N = 1;  RBmod_1  = model_project(HFmod,Vfull(:,1:N),Vfull(:,1:N));
N = 2;  RBmod_2  = model_project(HFmod,Vfull(:,1:N),Vfull(:,1:N));
N = 3;  RBmod_3  = model_project(HFmod,Vfull(:,1:N),Vfull(:,1:N));
N = 5;  RBmod_5  = model_project(HFmod,Vfull(:,1:N),Vfull(:,1:N));
N = 10; RBmod_10 = model_project(HFmod,Vfull(:,1:N),Vfull(:,1:N));    
%%
for i_fig = 1:2
    if i_fig == 1
        models = {ANNmod_g_1_12, ANNmod_3_12, ANNmod_5_16};
        tag = 'ANN';
    else    
        models = {RBmod_1, RBmod_3, RBmod_5, RBmod_10};
        tag = 'RB';
    end

    clear size_opt
    size_opt.subplot_width = 400;
    size_opt.subplot_heigth = 180;
    size_opt.margin_plot_v = 40;
    size_opt.margin_plot_h = 40;
    size_opt.bottom_margin = 40;
    size_opt.plot_legend = 1;
    win_sizes = plotting_get_axes(nrows,ncols,size_opt);
    fig_main = figure('units','pixel','position',[100 100 win_sizes.width win_sizes.heigth]);
    if do_inset
        win_sizes_mono = plotting_get_axes(1,1,struct('subplot_width',120,'subplot_heigth',80));
        for iMod = 1:length(ds)
            fig_mono{iMod} = figure('units','pixel','position',[100 100 win_sizes_mono.width win_sizes_mono.heigth]);
            subplot('Position',win_sizes_mono.axs{1}.coord_norm)
        end
    end
    cmap = get(0, 'DefaultAxesColorOrder');

    if do_plot_u
        leg = {'u_1','u_2','y (HF)'};
    else
        leg = {'y (HF)'};
    end

    for iMod = 1:length(models)
        leg = [leg {sprintf('y (%s, n = %d)',tag,models{iMod}.nX)}];
        ds_out = model_solve(ds,models{iMod});
        for k=1:1+do_inset        
            for i = 1:length(ds)
                if k == 1
                    figure(fig_main);
                    subplot('Position',win_sizes.axs{i}.coord_norm)
                else
                    figure(fig_mono{i});
                end
    %         subplot(nrows,ncols,i)
                if iMod == 1
                    if do_plot_u
                        plot(ds_out{i}.tt,ds_out{i}.uu,sty_u,'linewidth',lwidU); hold on;
                    end
                    plHF = plot(ds_out{i}.tt,ds_out{i}.yy_ex,sty_HF,'linewidth',lwidYHF); hold on;
                    cmap = get(0, 'DefaultAxesColorOrder');
                    if i == 1, legplots = plHF(1); end
                end
                plROM = plot(ds_out{i}.tt,ds_out{i}.yy,sty_ROM,'linewidth',lwidYANN,'color',cmap(iMod,:)); hold on;  
                if i == 1, legplots = [legplots plROM(1)]; end
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

    if do_inset
        for i = 1:length(ds)
            figure(fig_main);
            subplot('Position',win_sizes.axs{i}.coord_norm)
            in = inset{i};
        %     plot(in(1) + [0     in(3)],in(2) + [0     0    ], '-')
        %     plot(in(1) + [in(3) in(3)],in(2) + [0     in(4)], '-')
        %     plot(in(1) + [in(3) 0    ],in(2) + [in(4) in(4)], '-')
        %     plot(in(1) + [0     0    ],in(2) + [in(4) 0    ], '-')
            col = [1 1 1] * .3;
            plot([in(1) in(2)], [in(3) in(3)], '-', 'Color', col)
            plot([in(2) in(2)], [in(3) in(4)], '-', 'Color', col)
            plot([in(2) in(1)], [in(4) in(4)], '-', 'Color', col)
            plot([in(1) in(1)], [in(4) in(3)], '-', 'Color', col)
        end
    end

    % for i = 1:length(ds)
    %     %subplot(nrows,ncols,i)
    %     legend(leg)
    % end
    figure(fig_main);
    subplot('Position',win_sizes.axs{1}.coord_norm)
    myleg = legend(legplots,leg,'Orientation','horizontal');    
    set(myleg,'position', win_sizes.get_legend_coord_norm(myleg.Position), 'units', 'normalized');    

    pause(1e-16)
    if do_save
        create_directory_if_not_found('fig')
        figure(fig_main);
        print(sprintf('fig/TB_samples_test_main_%s',tag),'-depsc'); 
        if do_inset
            for i = 1:length(ds)
                figure(fig_mono{i});
                print(sprintf('fig/TB_samples_test_inset_%d',i),'-depsc'); 
            end
        end
    end
end