function [V,output] = POD_projection(X,opt)
    % Performs POD on the snapshots matrix X.
    % X: each column is a snapshot
    % V: each column is a basis
    
    fprintf('Computing POD through SVD decomposition...\n')
    
    if ~isfield(opt,'get_full_V')
        opt.get_full_V = 0;
    end
    if ~isfield(opt,'fixed_N')
        opt.fixed_N = 0;
    end
    if ~isfield(opt,'toll')
        opt.toll = 1e-2;
    end
    if ~isfield(opt,'do_plot')
        opt.do_plot = 0;
    end
    if ~isfield(opt,'show_epsilon')
        opt.show_epsilon = 0;
    end
        
    [~,S,V] = svd(X','econ');
    s = diag(S);
        
    if ~opt.get_full_V
        if opt.fixed_N
            N = opt.N;
        else
            sigma2_cum = cumsum(s.^2);
            [~,N] = max( opt.toll^2>=(1-sigma2_cum./sigma2_cum(end)) ) ;
        end

        eps = POD_getepsilon(s,N);

        fprintf('POD:  N   = %d\n',N)
        fprintf('      eps = %1.2e\n',eps)
        
        V = V(:,1:N);
        %U = U(:,1:N);    
    end
    
    if nargout > 1
        output.s = s;
        if ~opt.get_full_V
            output.N = N;
            output.eps = eps;
        end
    end
        
    if opt.do_plot
        figure()
        semilogy(s,'o-','linewidth',1);
%         if ~opt.get_full_V
%             hold on
%             semilogy(s(1:N),'o-','linewidth',2);
%         end
        hold on

        min_plot = min(s);
        max_plot = max(s);
        
        if opt.show_epsilon
            epss = [];
            for Ncurr = 1:length(s)
                epss = [epss POD_getepsilon(s,Ncurr)];
            end
            semilogy(epss,'o-','linewidth',1);
            legend('\sigma','\epsilon')

            min_plot = max(min(min_plot, min(epss)),1e-16);
            max_plot = max(max_plot, max(epss));
        end
        
        if ~opt.get_full_V
            semilogy([N N], [min_plot max_plot],'k--','linewidth',1);
            semilogy(s(1:N),'o-','linewidth',1);
        end
        
        grid on
    end
            
end