function model = testcase_2var_getmodels(problem, name)

    fprintf('building toy problem model...')
    model.problem = problem;
    

    model.nA = 1;
    model.alpha_min = 0;
    model.alpha_max = 1;
    
    switch name
        case 'exp'
            model.f_alpha = @(x,u,a) [a*x(1);0];
            model.dfdx = @(x,u,a) [a, 0;0, 0];
            model.dfda = @(x,u,a) [x(1);0];
    end

    model.nX = 2;
    model.x0_min = [1;1];
    model.x0_max = [1.1;3];
    model.dt = 1e-2;
    model.advance_type = 'nonlinear_explicit';
    model.u_implicit = 0;
    model.output_type = 'insidestate';
    
    fprintf(' done!\n')
    
end