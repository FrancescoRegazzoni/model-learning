clear
x0 = rand(1000,1);

%%
func.residuals = @(x) sin(x);
tic
x_sol = levenberg_marquardt(func,x0,struct('do_plot',0));
toc
norm(x_sol)
%%
func.residuals = @(x) sin(x);
func.jacobian = @(x) diag(cos(x));
tic
x_sol = levenberg_marquardt(func,x0,struct('do_plot',0));
toc
norm(x_sol)
%%
tic
x_sol = fsolve(@(x) sin(x),x0);
toc
norm(x_sol)