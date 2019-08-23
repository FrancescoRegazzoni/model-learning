clear
do_save = 0;

%% problem definition
problem = problem_get('thermalblock','TB9.ini');

%% Model generation
HFmod = problem.get_model(problem);

%% Load dataset
dataset_def.problem = problem;
dataset_def.type = 'file';
dataset_def.source = 'samples_x_rnd_A.mat;1:50';
ds = dataset_get(dataset_def);

%% Plots
ds_idx = [19 12 25 47];

fespace = HFmod.fespace;
n_vertices = size(fespace.mesh.vertices,1);
n1 = size(fespace.mesh.X,1);
n2 = size(fespace.mesh.X,2);
X = fespace.mesh.X;
Y = fespace.mesh.Y;
           
for iS = 1:length(ds_idx)
    figure('units','pixel','outerposition',[100 100 350 350]) 
    vec = ds{ds_idx(iS)}.xx(:,10);
    Z = reshape(vec(1:n_vertices),n1,n2);
    [C,h] = contourf(X,Y,Z,100);
    set(h,'LineColor','none')
    axis equal
    colorbar
    caxis([0 .45])
    hold on
    plot([1 3 5]*.25,[1 3 5]*.25,'o','MarkerEdgeColor','none', 'MarkerFaceColor',[.75 0 0])
     
    if do_save
        create_directory_if_not_found('fig')
        print(sprintf('fig/TB_snapshot_%d',iS),'-depsc'); 
    end
end

% figure('units','pixel','outerposition',[100 100 1000 300])            
% for iS = 1:length(ds_idx)
%     subplot(1,length(ds_idx),iS)
%     vec = ds{ds_idx(iS)}.xx(:,10);
%     Z = reshape(vec(1:n_vertices),n1,n2);
%     [C,h] = contourf(X,Y,Z,100);
%     set(h,'LineColor','none')
%     axis equal
% %     if iS == length(ds_idx)
%         colorbar
% %     end
% end
%%
%
% figure()
% imagesc(x,y,Z);
%
% figure()
% surf(fespace.mesh.X,fespace.mesh.Y,Z)
% hold on
% imagesc(fespace.mesh.X,fespace.mesh.Y,Z)

% %%
% 
% figure()
% plot_fe_function(ds{1}.xx(:,5),HFmod.fespace,'contourf',100)
% %set(axes,'LineColor','none')
% figure()
% plot_fe_function(ds{1}.xx(:,5),HFmod.fespace)
% 

