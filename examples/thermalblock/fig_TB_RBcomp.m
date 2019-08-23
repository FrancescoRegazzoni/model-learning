% NB: this script generates figures based on the previously computed
% comparison between ANNs and RB, stored in the file data_file. To
% regenerate the data file, run the script TB_RB_comparison.m.

clear

do_save = 0;
data_file = 'RB_comparison_data.mat';

sourceterm_pb = 0;
if sourceterm_pb
    problem = problem_get('thermalblock','TB9_sourceterm.ini');
else
    problem = problem_get('thermalblock','TB9.ini');
end
problem.goto_dir()
HFmod = problem.get_model(problem);
%%
RBev = load(data_file);
win_size = [500 500 300 400];
st_RB = '*-';
st_ANN = {'o-','+-'};

Nmax = max(RBev.optRBeval.N_values);
leg_vals = {'RB','ANN, out', 'ANN, in'};
figure('units','pixel','outerposition',win_size);
semilogy(RBev.optRBeval.N_values,RBev.ret.epss,'k--','linewidth',1.2); hold on
semilogy(RBev.optRBeval.N_values,RBev.ret.errs,st_RB,'linewidth',1.5); hold on
for iC = 1:2
    semilogy(RBev.ret.compmod{iC}.comp_nX,RBev.ret.compmod{iC}.comp_err,st_ANN{iC},'linewidth',1.5); hold on
end
% legend('\epsilon','error_{RB}','error_{ANN, ext}','error_{ANN, int}','location','southwest')
legend(['\epsilon_n',leg_vals]','location','southwest')
xlabel('n')
ylabel('test error')
grid on
if do_save, print('fig/TB_RBcomp_err','-depsc'); end

figure('units','pixel','outerposition',win_size);
semilogy(RBev.optRBeval.N_values,RBev.ret.tims/RBev.ret.timeRef,st_RB,'linewidth',1.5); hold on
for iC = 1:2
    semilogy(RBev.ret.compmod{iC}.comp_nX,RBev.ret.compmod{iC}.comp_tim/RBev.ret.timeRef,st_ANN{iC},'linewidth',1.5); hold on
end
% legend('time_{RB}/time_{HF}','time_{ANN, ext}/time_{HF}','time_{ANN, int}/time_{HF}','location','southeast')
legend(leg_vals,'location','southeast')
xlabel('n')
ylabel('time / time_{HF}')
grid on
axis([0 Nmax 3e-4 1e-2])
if do_save, print('fig/TB_RBcomp_tim','-depsc'); end

figure('units','pixel','outerposition',win_size);
loglog(RBev.ret.tims/RBev.ret.timeRef,RBev.ret.errs,st_RB,'linewidth',1.5); hold on
for iC = 1:2
    loglog(RBev.ret.compmod{iC}.comp_tim/RBev.ret.timeRef,RBev.ret.compmod{iC}.comp_err,st_ANN{iC},'linewidth',1.5); hold on
end
% legend('time_{RB}/time_{HF}','time_{ANN, ext}/time_{HF}','time_{ANN, int}/time_{HF}','location','southeast')
legend(leg_vals,'location','northwest')
xlabel('time / time_{HF}')
ylabel('test error')
grid on
axis([3e-4 1e-2 5e-3 5e-1])
if do_save
    create_directory_if_not_found('fig')
    print('fig/TB_RBcomp_timerr','-depsc'); 
end