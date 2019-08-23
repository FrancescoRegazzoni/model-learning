function D = numerical_diff_testconvergence(f,x0,opt)
    % Utility to check convergence of numerical differentiation of 'f'
    % around the point 'x0'.

    %% default values
    opt.dummy = 0;
    if ~isfield(opt,'idx_out')
        opt.idx_out = [];
    end
    if ~isfield(opt,'idx_in')
        opt.idx_in = [];
    end
    if ~isfield(opt,'jac_exact')
        opt.jac_exact = [];
    end
    
    ep_vals = 10.^(5:-1:-10);
    n_ep = length(ep_vals);
    
    %% execution
    n = size(x0,1);
    f0 = f(x0);
    m = size(f0,1);
    D = zeros(m,n,n_ep);
    
    if isempty(opt.idx_out)
        i_vals = 1:m;
    else
        i_vals = opt.idx_out;
    end
    ni = length(i_vals);
    if isempty(opt.idx_in)
        j_vals = 1:n;
    else
        j_vals = opt.idx_in;
    end
    nj = length(j_vals);
    myfig = figure('units','normalized','outerposition',[0 0 1 1]);
    for idx_j = 1:nj
        j = j_vals(idx_j);
        for e = 1:n_ep
            x = x0;
            ep = ep_vals(e);
            x(j) = x(j) + ep;
            D(:,j,e) = (f(x)-f0)/ep;
        end
        for idx_i = 1:ni
            i = i_vals(idx_i);
            subplot(ni,nj, nj*(idx_i-1) + idx_j)
            semilogx(ep_vals,squeeze(D(i,j,:)),'+-','linewidth',1.5)
            if ~isempty(opt.jac_exact)
                hold on
                plot([min(ep_vals) max(ep_vals)],opt.jac_exact(i,j)*[1 1],'r--','linewidth',1.2)
                set(gca,'xtick',[])
                set(gca,'ytick',[])
            end
        end
        pause(1e-16)
    end
        
end