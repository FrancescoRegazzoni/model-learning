function model = NTL1_getmodel_blackbox(problem)

    fprintf('building NTL1 model (blackbox version)...')
    
    % assign the problem struct to the model struct
    model.problem = problem;
    
    model.blackbox = 1;
    model.forward_function = @forward_function;

    fprintf(' done!\n')
    
    function output = forward_function(test, options)
        
        % constants
        N = 1e3;
        gamma = 40; 
        % integration time step
        dt = 5e-3; 
        % initial state
        x0 = zeros(N,1); 
        
        % initialization
        tt = test.tt(1):dt:test.tt(end);
        if isa(test.uu,'function_handle')
            ufunc = test.uu;
            uu = ufunc(tt);
        else
            uu = interp_time_series(test.tt,test.uu,tt);
        end
        nT = length(tt);
        yy = zeros(1,nT);
        
        % time loop
        for iT = 1:nT
            if iT==1
                x = x0;
            else
                dx = zeros(N,1);
                dx(1) = -2*x(1) + x(2) + 2 - exp(gamma*x(1)) - exp(gamma*(x(1)-x(2))) + uu(iT);
                dx(2:N-1) = -2*x(2:N-1) + x(1:N-2) + x(3:N) + exp(gamma*(x(1:N-2)-x(2:N-1))) - exp(gamma*(x(2:N-1)-x(3:N)));
                dx(N) = - x(N) + x(N-1) - 1 + exp(gamma*(x(N-1)-x(N)));
                x = x + dt*dx;
            end
            yy(1,iT) = x(1);
        end      
        
        % output
        output.tt = tt;
        output.uu = uu;
        output.yy = yy;
    end

end