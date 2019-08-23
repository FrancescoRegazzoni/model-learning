clear
% problem definition
problem = problem_get('transmission_line_circuit','INVC.ini');

% Model generation
HFmod = problem.get_model(problem);

% single test (see paper)
opt_solve.do_plot = 1;
figure();

test_solve.tt = [0 10];
test_solve.uu = @(t) .5*(1+cos(2*pi*t/10));
model_solve(test_solve,HFmod,opt_solve);


