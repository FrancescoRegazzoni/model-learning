clear
do_save = 0;

%%
problem = problem_get('thermalblock','TB9.ini');
problem.goto_dir()
HFmod = problem.get_model(problem);

%%
sourceterm_pb = 0;
load_ANNs;
%%
POD_dataset = 'samples_x_const_lhs10.mat;1:10|samples_x_rnd_A.mat;1:50';
X = build_snapshots_matrix(problem,POD_dataset);
%% get RB models
optPOD.get_full_V = 1;
[Vfull,outputPOD] = POD_projection(X,optPOD);
N = 1;  RBmod_1  = model_project(HFmod,Vfull(:,1:N),Vfull(:,1:N));
N = 2;  RBmod_2  = model_project(HFmod,Vfull(:,1:N),Vfull(:,1:N));
N = 3;  RBmod_3  = model_project(HFmod,Vfull(:,1:N),Vfull(:,1:N));
N = 10; RBmod_10 = model_project(HFmod,Vfull(:,1:N),Vfull(:,1:N));    

%% comparison plot
dataset_def.problem = problem;
dataset_def.type = 'file';
dataset_def.source = 'samples_rnd.mat;8';
test = dataset_get(dataset_def);

test_solve = test{1};

%outHF = model_solve(test_solve,HFmod);
outANN_1 = model_solve(test_solve,ANNmod_g_1_12);
outANN_2 = model_solve(test_solve,ANNmod_g_2_12);
outANN_3 = model_solve(test_solve,ANNmod_3_12);
outANN_4 = model_solve(test_solve,ANNmod_4_12);
outANN_5 = model_solve(test_solve,ANNmod_5_16);
outRB_1  = model_solve(test_solve,RBmod_1);
outRB_2  = model_solve(test_solve,RBmod_2);
outRB_3  = model_solve(test_solve,RBmod_3);
outRB_10 = model_solve(test_solve,RBmod_10);

%%
figure('units','pixel','outerposition',[500 500 700 400]);
cmap = get(0, 'DefaultAxesColorOrder');
lwid = 1.5;
plHF = plot(test_solve.tt,test_solve.yy,'k-','linewidth',lwid); hold on;
plR1 = plot(outANN_1.tt,outANN_1.yy,'--','linewidth',lwid,'color',cmap(1,:)); hold on;
plR2 = plot(outANN_2.tt,outANN_2.yy,'--','linewidth',lwid,'color',cmap(2,:)); hold on;
plR3 = plot(outANN_3.tt,outANN_3.yy,'--','linewidth',lwid,'color',cmap(3,:)); hold on;
plR4 = plot(outANN_3.tt,outANN_4.yy,'--','linewidth',lwid,'color',cmap(4,:)); hold on;
plR5 = plot(outANN_3.tt,outANN_5.yy,'--','linewidth',lwid,'color',cmap(5,:)); hold on;
xlabel('t')
ylabel('y')
legend([plHF(1) plR1(1) plR2(1) plR3(1) plR4(1) plR5(1)],{'HF model (N=3721)','ANN model (n=1, ext)','ANN model (n=2, ext)','ANN model (n=3, int)','ANN model (n=4, int)','ANN model (n=5, int)'},'location','eastoutside')
% legend([plHF,plR1,plR2,plR3],{'HF model (N=3721)','ANN model (n=1)','ANN model (n=2)','ANN model (n=3)'},'location','north')
if do_save
    print('TB_comparison_ANN','-depsc'); 
end
%%
figure('units','pixel','outerposition',[500 500 700 400]);
cmap = get(0, 'DefaultAxesColorOrder');
lwid = 1.5;
plHF = plot(test_solve.tt,test_solve.yy,'k-','linewidth',lwid); hold on;
plR4 = plot(outRB_1.tt ,outRB_1.yy ,'--','linewidth',lwid,'color',cmap(1,:)); hold on;
plR5 = plot(outRB_2.tt ,outRB_2.yy ,'--','linewidth',lwid,'color',cmap(2,:)); hold on;
plR6 = plot(outRB_3.tt ,outRB_3.yy ,'--','linewidth',lwid,'color',cmap(3,:)); hold on;
plR7 = plot(outRB_10.tt,outRB_10.yy,'--','linewidth',lwid,'color',cmap(5,:)); hold on;
xlabel('t')
ylabel('y')
legend([plHF(1) plR4(1) plR5(1) plR6(1) plR7(1)],{'HF model (N=3721)','RB model (n=1)','RB model (n=2)','RB model (n=3)','RB model (n=10)'},'location','eastoutside')
if do_save
    print('TB_comparison_RB','-depsc'); 
end