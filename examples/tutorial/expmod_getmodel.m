function model = expmod_getmodel(problem)

    fprintf('building toy problem model...')
    
    model.problem = problem;
    
    model.nX = 1; % number of internal states
    model.x0_min = 1;
    model.x0_max = 1.1;
    model.nA = 1; % number of parameters
    model.alpha_min = 0;
    model.alpha_max = 1;
    model.dt = 1e-2; % integration time step
    
    % Model dynamics
    model.advance_type = 'nonlinear_explicit';
    model.f_alpha = @(x,u,a) a*x;
    model.dfdx = @(x,u,a) a;
    model.dfda = @(x,u,a) x;
    
    % output definition
    model.output_type = 'insidestate';
    
    fprintf(' done!\n')
    
end