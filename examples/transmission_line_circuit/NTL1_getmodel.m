function model = NTL1_getmodel(problem,N)

    fprintf('building NTL1 model...')
    
    model.problem = problem;
    
    gamma = 40;
    
    G = sparse(1,N);
    G(1,1) = 1;
    
    model.nX = N;
    model.x0 = zeros(N,1);
    model.dt = 5e-3;
    model.advance_type = 'nonlinear_explicit';
    model.u_implicit = 0;
    model.output_type = 'linear';
    model.G = G;
    model.g0 = 0;
    model.f = @f;  
    
    fprintf(' done!\n')
    
    function dx = f(x,u)
        dx = zeros(N,1);
        dx(1) = -2*x(1) + x(2) + 2 - exp(gamma*x(1)) - exp(gamma*(x(1)-x(2))) + u;
        dx(2:N-1) = -2*x(2:N-1) + x(1:N-2) + x(3:N) + exp(gamma*(x(1:N-2)-x(2:N-1))) - exp(gamma*(x(2:N-1)-x(3:N)));
        dx(N) = - x(N) + x(N-1) - 1 + exp(gamma*(x(N-1)-x(N)));
    end

end