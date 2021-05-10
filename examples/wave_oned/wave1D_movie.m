clear

do_save = 1;
%%
problem = problem_get('wave_oned','wave1D.ini');
problem.goto_dir()
HFmod = problem.get_model(problem);

%% test
% test.tt = [0 10];
% test.uu = @(t) [sin(2*pi/.1*t); cos(2*pi/.2*t)];

% dataset_def.problem = problem;
% dataset_def.type = 'file';
% dataset_def.source = 'samples_x_rnd.mat;1';
% ds = dataset_get(dataset_def);
% test = ds{1};

test.tt = 0:HFmod.dt:5;
optGen.dim = problem.nU;
optGen.umin = problem.u_min;
optGen.umax = problem.u_max;
optGen.method = 'coloured_noise'; 
optGen.time_scale = .05;
rng(1)
test.uu = get_random_time_course(test.tt,optGen);

figure()
out = model_solve(test,HFmod,struct('save_x',1,'do_plot',1));
%%
figure()
n = size(out.xx,1)/2;
for iT = round(linspace(1,length(out.tt),10))
    subplot(2,1,1)
    plot(out.xx(1:n+1,iT),'-','linewidth',2)
    hold on
    subplot(2,1,2)
    plot(out.xx(n+2:end,iT),'-','linewidth',2)
    hold on
end
%%
if do_save
    create_directory_if_not_found('fig');
    create_directory_if_not_found('fig/wave1D_video');
end

moviefig = figure('units','pixel','outerposition',[100 100 350 300]);
n = size(out.xx,1)/2;
for iT = 2:length(out.tt)
    %figure(moviefig)
    plot(out.xx(1:n+1,iT),'-','linewidth',2)
    hold on
    plot([1;n+1],out.uu(:,iT),'o','MarkerSize',7,'MarkerEdgeColor','none', 'MarkerFaceColor',[0 .75 0])
    plot(n/2,out.xx(n/2,iT),'o','MarkerSize',7,'MarkerEdgeColor','none', 'MarkerFaceColor',[.75 0 0])
    axis([-inf inf problem.y_min,problem.y_max])
    title(sprintf('t = %1.2f',out.tt(iT)));
    %set(gca,'Visible','off')
%     color = get(hFig,'Color');
%     set(gca,'XColor',color,'YColor',color,'TickDir','out')
%     set(gca,'box','off');
    set(gca,'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[])
    hold off
    pause(1e-16)
    if do_save
        print(sprintf('fig/wave1D_video/wave1D_frame_%06d',iT),'-dpng'); 
    end
end
