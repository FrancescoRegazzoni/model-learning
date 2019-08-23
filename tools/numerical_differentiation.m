function D = numerical_differentiation(f,x0,opt)
    % Approximates the first derivative of 'f' around 'x0'.

    %% default values
    opt.dummy = 0;
    if ~isfield(opt,'ep_rel')
        opt.ep_rel = 1e-6;
    end
    if ~isfield(opt,'ep_min')
        opt.ep_min = 1e-10;
    end
    opt.ep_min = 1e-10;
    
    %% execution
    n = size(x0,1);
    f0 = f(x0);
    D = zeros(size(f0,1),n);
    for i = 1:n
        x = x0;
        ep = max(opt.ep_min,abs(x(i))*opt.ep_rel);
        x(i) = x(i) + ep;
        D(:,i) = (f(x)-f0)/ep;		
    end
        
end