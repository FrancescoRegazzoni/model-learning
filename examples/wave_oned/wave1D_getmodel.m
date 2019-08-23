function model = TB_getmodel(problem, input_type)

    fprintf('building 1-dimensional wave equation model ...')
    
    model.problem = problem;
    
    L = 1;
    h = 1e-2;
    dt = 1e-2;
    c = 2;
    
    n = round(L/h);
    if abs(n*h-L)/abs(L) > 1e-12
        error('L must be multiple of h')
    end
    
    stiff_matrix = c^2*(2*eye(n+1)-diag(ones(n,1),1)-diag(ones(n,1),-1))/h^2;
        
%     model.Kv0 = [eye(n+1)/dt     zeros(n+1,n-1); ...
%                  zeros(n-1,n+1)  eye(n-1)/dt   ];

    Z = zeros(2*n);
    z = zeros(2*n,1);
    
    model.Kf0 = [zeros(n+1)               [zeros(1,n-1); -eye(n-1); zeros(1,n-1)]; ...
                 stiff_matrix(2:end-1,:)  zeros(n-1)                            ];
    model.Kv0 = eye(2*n);
    model.Bf0 = Z;
    model.Bv0 = eye(2*n);
    model.b0 = z;    
    for iU = 1:model.problem.nU         
        model.Kf{iU} = Z;  
        model.Kv{iU} = Z;
        model.Bf{iU} = Z;
        model.Bv{iU} = Z;
        model.b{iU}  = z;
    end
    
    %BC
    model.Kf0(1,1) = 1;
    model.Kv0(1,1) = 0;
    model.Bv0(1,1) = 0;
    model.b{1}(1) = 1;
    model.Kf0(n+1,n+1) = 1;
    model.Kv0(n+1,n+1) = 0;
    model.Bv0(n+1,n+1) = 0;
    model.b{2}(n+1) = 1;
    
    model.Kv0 = sparse(model.Kv0);
    model.Kf0 = sparse(model.Kf0);
    model.Bv0 = sparse(model.Bv0);
    model.Bf0 = sparse(model.Bf0);
    for iU = 1:model.problem.nU 
        model.Kv{iU}  = sparse(model.Kv{iU} );
        model.Kf{iU}  = sparse(model.Kf{iU} );
        model.Bv{iU}  = sparse(model.Bv{iU} );
        model.Bf{iU}  = sparse(model.Bf{iU} );
    end
    
    g0 = zeros(model.problem.nY,1);
    G = sparse(model.problem.nY,2*n); 
    G(1,round(n/2+1)) = 1;
    
    model.nX = 2*n;
    model.x0 = zeros(2*n,1);
    model.dt = dt;
    model.advance_type = 'linear_advance_timeexplicit_uaffine';
    model.u_implicit = 1;
    model.output_type = 'linear';
    model.G = G;
    model.g0 = g0;
    model.make_plot = @make_plot;
            
    fprintf('done!\n')
    
    function make_plot(x,opt)
        
        opt.dummy = 0;
        if ~isfield(opt,'new_figure')
            opt.new_figure = 0;
        end
        
        if opt.new_figure
            figure('units','normalized','outerposition',[0 0 .3 .6]);
        end

        subplot(2,1,1)
        plot(x(1:n+1),'o-','linewidth',2)
%         axis([-inf inf problem.y_min,problem.y_max])
        axis([-inf inf -inf inf])
        subplot(2,1,2)
        plot(x(n+2:end),'o-','linewidth',2)
        axis([-inf inf -inf inf])    
    end

end