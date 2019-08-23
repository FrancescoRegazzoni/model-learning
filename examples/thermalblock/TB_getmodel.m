function model = TB_getmodel(problem, input_type)

    % - \grad \cdot (mu(x) \cdot u) = f(x)
    % u = 0
    % u = neumann_val
    
    if nargin == 1
        input_type = 0; % 0 = mu, 1 = f
    end
    fprintf('building Thermal Block (TB) model:\n')
    
    if ~exist('create_mesh','file')
        error('feamat library not found. Please read the ''README'' file.')
    end
    
    model.problem = problem;
    
    muMin = model.problem.u_min(1);
    muMax = model.problem.u_max(1);
    if input_type == 0
        neumann_val = 10;
        full_dirichlet = 0;
    else
        neumann_val = 0;
        full_dirichlet = 0;
    end
    param_mu_000 = 10;
    mu_constant = 10; 
    param_f_000 = 0;
    probesX = [ .25 .75 1.25 ];
    probesY = [ .25 .75 1.25 ];

    % geometry  definition
    bottom_left_corner_x = 0;
    bottom_left_corner_y = 0;

    L = 1.5;
    H = 1.5;

    % number of elements
    n_elements_x = 30;
    n_elements_y = 30;

    mesh = create_mesh(bottom_left_corner_x, ...
                       bottom_left_corner_y, ...
                       L,H,n_elements_x,n_elements_y);

    % boundary conditions   
    if full_dirichlet
        bc_flags = [0 0 0 0];
    else
        bc_flags = [0 0 1 0];
    end
    dirichlet_functions = @(x) [0;0;0;0];
    neumann_functions = @(x) [neumann_val;0;0;0];

    % finite elements space definition
    fespace = create_fespace(mesh,'P2',bc_flags);

    % forcing term
    f = @(x) 0*x(1,:);

    % thermal block parameters
    %param = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 1 ];


    n_vertices = size(fespace.mesh.vertices,1);
    n1 = size(fespace.mesh.X,1);
    n2 = size(fespace.mesh.X,2);


    nProbes = length(probesX);
    for iProbe = 1:nProbes
        [~,idxProbeX(iProbe)]=min(abs(fespace.mesh.X(:,1)-probesX(iProbe)));
        [~,idxProbeY(iProbe)]=min(abs(fespace.mesh.Y(1,:)-probesY(iProbe)));
    end

    fprintf('   assembling matrices...\n')
    M = assemble_mass(fespace);
    nFEM = size(M,1);    
    Z = zeros(nFEM,nFEM);
    D = apply_dirichlet_bc_matrix(Z,fespace,1);
    model.Kv0 = apply_dirichlet_bc_matrix(M,fespace,1) - D;
    if input_type == 0
        model.Kf0 = D;
    else
        model.Kf0 = apply_dirichlet_bc_matrix(assemble_stiffness(mu_constant,fespace),fespace,1);
    end
    model.Bv0 = model.Kv0;
    model.Bf0 = Z;
    model.b0  = apply_dirichlet_bc_rhs( ...
                   apply_neumann_bc( ...
                      assemble_rhs( ...
                         fespace,f) ...
                      ,fespace,neumann_functions) ...
                   ,fespace,dirichlet_functions);

    for iU = 1:model.problem.nU 
        fprintf('   assembling matrices (parameter %d of %d)...\n',iU,model.problem.nU )
        u = zeros(model.problem.nU,1);
        u(iU) = 1;
        if input_type == 0
            mu = get_piecewise_function(u,param_mu_000);
            model.Kf{iU} = apply_dirichlet_bc_matrix(assemble_stiffness(mu,fespace),fespace,1) - D;
            model.b{iU}  = zeros(nFEM,1);
        else
            f = get_piecewise_function(u,param_f_000);
            model.Kf{iU} = Z;
            model.b{iU} = apply_dirichlet_bc_rhs(assemble_rhs(fespace,f),fespace,dirichlet_functions); % being homogeneous I do not need to set to zero the dirichlet BC
        end
        model.Kv{iU} = Z;
        model.Bf{iU} = Z;
        model.Bv{iU} = Z;
    end

    fprintf('   making sparse the matrices...\n')
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
    
    fprintf('   building the output matrix...\n')
    G = sparse(model.problem.nY,nFEM); 
    for iProbe = 1:nProbes
        G(iProbe,(idxProbeY(iProbe)-1)*n1+idxProbeX(iProbe)) = 1;
    end
    
    model.fespace = fespace;
        
    model.nX = nFEM;
%     model.nU = nParams;
%     model.nY = 3;
    model.x0 = zeros(nFEM,1);
%     model.umin = ones(model.nU,1)*muMin;
%     model.umax = ones(model.nU,1)*muMax;
    model.dt = 1e-2;
%     model.T = 1;
%     model.advance_type = 'solveonestep_linear';
%     model.solveonestep_linear = @solveonestep_linear;
    model.advance_type = 'linear_advance_timeexplicit_uaffine';
    model.u_implicit = 1;
%     model.output_type = 'nonlinear';
%     model.g = @g;   
    model.output_type = 'linear';
    model.G = G;
    model.g0 = zeros(model.problem.nY,1);
    model.make_plot = @make_plot;
            
    fprintf('done!\n')
    
%     function [K,rhs] = solveonestep_linear(x,u,dt)
%         % returns K and rhs such that:
%         %       K(u,dt) x^k+1 = rhs(x^k,u,dt)
%         
%         if nargin == 2
%             dt = model.dt;
%         end
%         
%         mu = get_mu(u);
%         A = assemble_stiffness(mu,fespace);
%         K = (M/dt)+A;
%         rhs = (M/dt)*x;    
%         rhs = rhs + assemble_rhs(fespace,f);
%         rhs = apply_neumann_bc(rhs,fespace,neumann_functions); 
%         [K,rhs] = apply_dirichlet_bc(K,rhs,fespace,dirichlet_functions);        
%     end

%     function xnew = solveonestep(x,u,dt)
%         [K,rhs] = solveonestep_linear(x,u,dt);
%         xnew = K \ rhs;
% %         if nargin == 2
% %             dt = mod.dt;
% %         end
% %         
% %         mu = get_mu(u);
% %         A = assemble_stiffness(mu,fespace);
% %         K = (M/dt)+A;
% %         rhs = (M/dt)*x;    
% %         rhs = rhs + assemble_rhs(fespace,f);
% %         rhs = apply_neumann_bc(rhs,fespace,neumann_functions); 
% %         [K,rhs] = apply_dirichlet_bc(K,rhs,fespace,dirichlet_functions);
% %         xnew = K \ rhs;
%     end

%     function y = g(x)    
%         solMat = reshape(x(1:n_vertices),n1,n2);
%         for iP = 1:nProbes
%             y(iP,1) = solMat(idxProbeX(iP),idxProbeY(iP));
%         end
%     end

    function ret = get_piecewise_function(param,base_param)    
        if length(param) == 9
            ret = @(x) param(1)*(x(1,:)<0.5).*(x(2,:)<0.5) ...
               + param(2)*(x(1,:)>=0.5).*(x(1,:)<1.0).*(x(2,:)>=0.0).*(x(2,:)<0.5) ...
               + param(3)*(x(1,:)>=1.0).*(x(1,:)<=1.5).*(x(2,:)>=0.0).*(x(2,:)<0.5) ...
               + param(4)*(x(1,:)>=0.0).*(x(1,:)<0.5).*(x(2,:)>=0.5).*(x(2,:)<1.0) ...
               + param(5)*(x(1,:)>=0.5).*(x(1,:)<1.0).*(x(2,:)>=0.5).*(x(2,:)<1.0) ...
               + param(6)*(x(1,:)>=1.0).*(x(1,:)<=1.5).*(x(2,:)>=0.5).*(x(2,:)<1.0) ...
               + param(7)*(x(1,:)>=0.0).*(x(1,:)<0.5).*(x(2,:)>=1.0).*(x(2,:)<=1.5) ...
               + param(8)*(x(1,:)>=0.5).*(x(1,:)<1.0).*(x(2,:)>=1.0).*(x(2,:)<=1.5) ...
               + param(9)*(x(1,:)>=1.0).*(x(1,:)<=1.5).*(x(2,:)>=1.0).*(x(2,:)<=1.5) ;
        elseif length(param) == 3
            ret = @(x) param(1)*(x(1,:)<0.5).*(x(2,:)<0.5) ...
               + base_param*(x(1,:)>=0.5).*(x(1,:)<1.0).*(x(2,:)>=0.0).*(x(2,:)<0.5) ...
               + param(2)*(x(1,:)>=1.0).*(x(1,:)<=1.5).*(x(2,:)>=0.0).*(x(2,:)<0.5) ...
               + base_param*(x(1,:)>=0.0).*(x(1,:)<0.5).*(x(2,:)>=0.5).*(x(2,:)<1.0) ...
               + param(3)*(x(1,:)>=0.5).*(x(1,:)<1.0).*(x(2,:)>=0.5).*(x(2,:)<1.0) ...
               + base_param*(x(1,:)>=1.0).*(x(1,:)<=1.5).*(x(2,:)>=0.5).*(x(2,:)<1.0) ...
               + param(2)*(x(1,:)>=0.0).*(x(1,:)<0.5).*(x(2,:)>=1.0).*(x(2,:)<=1.5) ...
               + base_param*(x(1,:)>=0.5).*(x(1,:)<1.0).*(x(2,:)>=1.0).*(x(2,:)<=1.5) ...
               + param(1)*(x(1,:)>=1.0).*(x(1,:)<=1.5).*(x(2,:)>=1.0).*(x(2,:)<=1.5) ;
        end
    end

    function mu = get_mu(param)        
        if length(param) == 9
            mu = @(x) param(1)*(x(1,:)<0.5).*(x(2,:)<0.5) ...
               + param(2)*(x(1,:)>=0.5).*(x(1,:)<1.0).*(x(2,:)>=0.0).*(x(2,:)<0.5) ...
               + param(3)*(x(1,:)>=1.0).*(x(1,:)<=1.5).*(x(2,:)>=0.0).*(x(2,:)<0.5) ...
               + param(4)*(x(1,:)>=0.0).*(x(1,:)<0.5).*(x(2,:)>=0.5).*(x(2,:)<1.0) ...
               + param(5)*(x(1,:)>=0.5).*(x(1,:)<1.0).*(x(2,:)>=0.5).*(x(2,:)<1.0) ...
               + param(6)*(x(1,:)>=1.0).*(x(1,:)<=1.5).*(x(2,:)>=0.5).*(x(2,:)<1.0) ...
               + param(7)*(x(1,:)>=0.0).*(x(1,:)<0.5).*(x(2,:)>=1.0).*(x(2,:)<=1.5) ...
               + param(8)*(x(1,:)>=0.5).*(x(1,:)<1.0).*(x(2,:)>=1.0).*(x(2,:)<=1.5) ...
               + param(9)*(x(1,:)>=1.0).*(x(1,:)<=1.5).*(x(2,:)>=1.0).*(x(2,:)<=1.5) ;
        elseif length(param) == 3
            mu = @(x) param(1)*(x(1,:)<0.5).*(x(2,:)<0.5) ...
               + param_mu_000*(x(1,:)>=0.5).*(x(1,:)<1.0).*(x(2,:)>=0.0).*(x(2,:)<0.5) ...
               + param(2)*(x(1,:)>=1.0).*(x(1,:)<=1.5).*(x(2,:)>=0.0).*(x(2,:)<0.5) ...
               + param_mu_000*(x(1,:)>=0.0).*(x(1,:)<0.5).*(x(2,:)>=0.5).*(x(2,:)<1.0) ...
               + param(3)*(x(1,:)>=0.5).*(x(1,:)<1.0).*(x(2,:)>=0.5).*(x(2,:)<1.0) ...
               + param_mu_000*(x(1,:)>=1.0).*(x(1,:)<=1.5).*(x(2,:)>=0.5).*(x(2,:)<1.0) ...
               + param(2)*(x(1,:)>=0.0).*(x(1,:)<0.5).*(x(2,:)>=1.0).*(x(2,:)<=1.5) ...
               + param_mu_000*(x(1,:)>=0.5).*(x(1,:)<1.0).*(x(2,:)>=1.0).*(x(2,:)<=1.5) ...
               + param(1)*(x(1,:)>=1.0).*(x(1,:)<=1.5).*(x(2,:)>=1.0).*(x(2,:)<=1.5) ;
        end
    end

    function make_plot(vecsol,opt)
        
        opt.dummy = 0;
        if ~isfield(opt,'new_figure')
            opt.new_figure = 0;
        end
        
        if opt.new_figure
            figure('units','normalized','outerposition',[0 0 .3 .6]);
        end
        
        plot_fe_function(vecsol,fespace)
        axis( [0 1.5 0 1.5 0 .3] )
    end

end