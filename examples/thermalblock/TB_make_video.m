clear
do_save = 1;

%% problem definition
problem = problem_get('thermalblock','TB9.ini');

%% Model generation
HFmod = problem.get_model(problem);

%% Load dataset
dataset_def.problem = problem;
dataset_def.type = 'file';
dataset_def.source = 'samples_x_detailed_rnd.mat;:';
ds = dataset_get(dataset_def);

%% Plots
ds_idx = 1;

fespace = HFmod.fespace;
n_vertices = size(fespace.mesh.vertices,1);
n1 = size(fespace.mesh.X,1);
n2 = size(fespace.mesh.X,2);
X = fespace.mesh.X;
Y = fespace.mesh.Y;

figure('units','pixel','outerposition',[100 100 350 350])
disp('press enter to go...')
pause()
for iT = 1:length(ds{ds_idx}.tt) 
    vec = ds{ds_idx}.xx(:,iT);
    Z = reshape(vec(1:n_vertices),n1,n2);
    [C,h] = contourf(X,Y,Z,100);
    set(h,'LineColor','none')
    axis equal
    colorbar
    caxis([0 .45])
    hold on
    plot([1 3 5]*.25,[1 3 5]*.25,'o','MarkerEdgeColor','none', 'MarkerFaceColor',[.75 0 0])
    hold off
    title(sprintf('t = %1.2f',ds{ds_idx}.tt(iT)));
    pause(1e-16)

    if do_save
        create_directory_if_not_found('fig');
        create_directory_if_not_found('fig/TB_video');
        print(sprintf('fig/TB_video/TB_frame_%06d',iT),'-dpng'); 
    end
end

%%
make_video('fig/TB_video', 'fig/TB_video.avi');

