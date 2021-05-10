function model = param_pend_getmodels(problem, type)

    fprintf('building parametric pendulum model...')
    model.problem = problem;
    
    c = 1;
    beta = 5;
    
    switch type
        case 'th0'
            model.nA = 1;
            model.alpha_min = -.2;
            model.alpha_max =  .2;
            model.f_alpha = @(x,u,a) [x(2);-c*x(2)-beta*sin(x(1)-a)];
            model.dfdx = @(x,u,a) [0                , 1; ...
                                   -beta*cos(x(1)-a), -c];
            model.dfda = @(x,u,a) [0; beta*cos(x(1)-a)];
        case 'th0_c'
            model.nA = 2;
            model.alpha_min = [-.2; .5];
            model.alpha_max = [.2; 1.5];
            model.f_alpha = @(x,u,a) [x(2);-a(2)*x(2)-beta*sin(x(1)-a(1))];
            model.dfdx = @(x,u,a) [0                , 1; ...
                                   -beta*cos(x(1)-a(1)), -a(2)];
            model.dfda = @(x,u,a) [0 , 0; ...
                                   beta*cos(x(1)-a(1)), -x(2)];
        otherwise
            error('non existent type')
    end

    model.nX = 2;
    model.x0_min = [-1;-.1];
    model.x0_max = [1; .1];
    model.dt = 1e-2;
    model.advance_type = 'nonlinear_explicit';
    model.u_implicit = 0;
    model.output_type = 'insidestate';
    
    fprintf(' done!\n')
    
end