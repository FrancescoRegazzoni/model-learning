function model = testcase_1var_getmodels(problem, name, useCN)

    fprintf('building toy problem model...')
    model.problem = problem;
    
    if nargin < 3
        useCN = 0;
    end

    if isequal(name, 'expstab_2a')
        model.nA = 2;
        model.alpha_min = [0;1];
        model.alpha_max = [1;3];
    else
        model.nA = 1;
        model.alpha_min = 0;
        model.alpha_max = 1;
    end
    switch name
        case 'linear'
            if useCN
                model.A_alpha = @(a) 0;
                model.b0_alpha = @(a) 2*a;
                model.dA_dalpha = @(a) 0;
                model.db0_dalpha = @(a) 2;
            else
                model.f_alpha = @(x,u,a) 2*a;
                model.dfdx = @(x,u,a) 0;
                model.dfda = @(x,u,a) 2;
            end
        case 'exp'
            if useCN
                model.A_alpha = @(a) a;
                model.b0_alpha = @(a) 0;
                model.dA_dalpha = @(a) 1;
                model.db0_dalpha = @(a) 0;
            else
                model.f_alpha = @(x,u,a) a*x;
                model.dfdx = @(x,u,a) a;
                model.dfda = @(x,u,a) x;
            end
        case 'expstab'
            if useCN
                model.A_alpha = @(a) -2*a;
                model.b0_alpha = @(a) 6*a;
                model.dA_dalpha = @(a) -2;
                model.db0_dalpha = @(a) 6;
            else
                model.f_alpha = @(x,u,a) 2*a*(3 - x);
                model.dfdx = @(x,u,a) -2*a;
                model.dfda = @(x,u,a) 2*(3 - x);
            end
        case 'bistab'
            if useCN
                error('CN not available')
            else
                model.f_alpha = @(x,u,a) 2*a*x*(3 - x);
                model.dfdx = @(x,u,a) 6*a - a*x;
                model.dfda = @(x,u,a) 2*x*(3 - x);
            end
        case 'expstab_2a'
            if useCN
                error('CN not available')
            else
                model.f_alpha = @(x,u,a) 2*a(1)*x*(a(2) - x);
                model.dfdx = @(x,u,a) 2*a(1)*(a(2) - 2*x);
                model.dfda = @(x,u,a) [2*x*(a(2) - x), 2*a(1)*x];
            end
    end

    model.nX = 1;
    model.x0_min = 1;
    model.x0_max = 1.1;
    model.dt = 1e-2;
    if useCN
        model.advance_type = 'linear_CN_uaffine';
    else
        model.advance_type = 'nonlinear_explicit';
    end
    model.u_implicit = 0;
    model.output_type = 'insidestate';
    
    fprintf(' done!\n')
    
end