function model = INVC_getmodel(problem,N)

    fprintf('building INVC model...')
    
    model.problem = problem;
    
    Vdd = 1;
    A = 5;
    f_inv = @(x) Vdd*tanh(-A*x);
    
    G = sparse(1,N);
    G(1,N) = 1;
    
    model.nX = N;
    model.x0 = zeros(N,1);
    model.dt = 1e-1;
    model.advance_type = 'nonlinear_explicit';
    model.u_implicit = 0;
    model.output_type = 'linear';
    model.G = G;
    model.g0 = 0;
    model.f = @f;  
    
    fprintf(' done!\n')
    
    function dx = f(x,u)
        dx = zeros(N,1);
        dx(1) = -x(1) + u;
        dx(2:N) = -x(2:N) + f_inv(x(1:N-1));
    end

end