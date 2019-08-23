function model = pendulum_getmodel(problem)

    fprintf('building pendulum model...')
    
    model.problem = problem;
    
    m = 1;
    c = 3;
    P = 2;
    L = 1;
    dt = 1e-2;
        
    model.nX = 2;
    model.x0 = zeros(2,1);
    model.dt = dt;
    
    model.advance_type = 'nonlinear_explicit';
    model.u_implicit = 0;
    model.f = @(x,u) [x(2); (u(2)*cos(x(1)) + (u(1)-P)*sin(x(1)))/m - c*x(2)];
    
    model.output_type = 'nonlinear';
    model.g = @(x) L*sin(x(1));
    
    fprintf(' done!\n')

end