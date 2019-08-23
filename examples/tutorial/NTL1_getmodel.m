function model = NTL1_getmodel(problem)

    fprintf('building NTL1 model...')
    
    % assign the problem struct to the model struct
    model.problem = problem;
    
    % constants
    N = 1e3;
    gamma = 40;   
    
    model.nX = N; % number of internal states
    model.x0 = zeros(N,1); % initial state
    model.dt = 5e-3; % integration time step
        
    % definition of model dynamics
    model.advance_type = 'nonlinear_explicit';
    model.f = @f; % defined below
    
    % definition of model output
    model.output_type = 'insidestate';
    
    fprintf(' done!\n')
    
    function dx = f(x,u)
        dx = zeros(N,1);
        dx(1) = -2*x(1) + x(2) + 2 - exp(gamma*x(1)) - exp(gamma*(x(1)-x(2))) + u;
        dx(2:N-1) = -2*x(2:N-1) + x(1:N-2) + x(3:N) + exp(gamma*(x(1:N-2)-x(2:N-1))) - exp(gamma*(x(2:N-1)-x(3:N)));
        dx(N) = - x(N) + x(N-1) - 1 + exp(gamma*(x(N-1)-x(N)));
    end

end