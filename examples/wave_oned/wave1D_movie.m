clear

problem = problem_get('wave_oned','wave1D.ini');
problem.goto_dir()
HFmod = problem.get_model(problem);

%% test
test.tt = [0 10];
test.uu = @(t) [sin(2*pi/.1*t); cos(2*pi/.2*t)];
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
moviefig = figure();
n = size(out.xx,1)/2;
for iT = round(linspace(1,length(out.tt),5000))
    figure(moviefig)
    plot(out.xx(1:n+1,iT),'-','linewidth',2)
    axis([-inf inf problem.y_min,problem.y_max])
    pause(1e-16)
end
