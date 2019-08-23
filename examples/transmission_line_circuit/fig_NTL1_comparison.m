clear
do_save = 0;
%% load problem and HF model
problem = problem_get('transmission_line_circuit','NTL1.ini');
problem.goto_dir()
HFmod = problem.get_model(problem);
%% load reduced models
load_ANNs;
%% Comparison plot
test_solve.tt = [0 1];
test_solve.uu = @(t) .5*(1+cos(2*pi*t));

outHF = model_solve(test_solve,HFmod);
outANN_1 = model_solve(test_solve,ANNmod_1_8_g);
outANN_2 = model_solve(test_solve,ANNmod_2_8_g);
outANN_3 = model_solve(test_solve,ANNmod_3_8_g);

figure('units','pixel','outerposition',[200 200 600 400]);
lwid = 1.5;
plot(outHF.tt,outHF.yy,'k','linewidth',lwid); hold on;
plot(outANN_1.tt,outANN_1.yy,'--','linewidth',lwid); hold on;
plot(outANN_2.tt,outANN_2.yy,'--','linewidth',lwid); hold on;
plot(outANN_3.tt,outANN_3.yy,'--','linewidth',lwid); hold on;
xlabel('t')
ylabel('y')
legend('HF model (N=1000)','ANN model (n=1)','ANN model (n=2)','ANN model (n=3)','location','north')
if do_save
    create_directory_if_not_found('fig')
    print('fig/NTL1_comparison','-depsc');
end
