function [y,x] = model_get_steady_state(model,u_steady,opt)
    
    % NB: this function is generated automatically by the script
    % "build_model_get_steady_state.m", for efficiency reasons. Do not 
    % modify it directly. To modify it:
    %  - enter into folder model_solve_builder;
    %  - modifiy the desired section (model_get_steady_state_header, 
    %    model_get_steady_state_loop_footer,etc.);
    %  - run "build_model_get_steady_state.m"
    %  - copy the script "model_get_steady_state.m" in this folder into "core"

    %% options
    opt.dummy = 0;
    if ~isfield(opt,'eps_x')
        opt.eps_x = 1e-4;
    end
    if ~isfield(opt,'eps_y')
        opt.eps_y = 1e-4;
    end
    if ~isfield(opt,'x0')
        opt.x0 = [];
    end
    
    %% initialization
    if isempty(opt.x0)
        x = model.x0;
    else
        x = opt.x0;
    end
    
    dt = model.dt;
    
    switch model.output_type
        case 'linear' 
            g = @(x) model.G*x+model.g0;
        case 'nonlinear'
            g = model.g;
        case 'insidestate'
            g = @(x) x(1:model.problem.nY);
        otherwise
            error('unknown output type %s',model.output_type)
    end
    
    y = g(x);
    
    y_norm = model.problem.y_max-model.problem.y_min;
    
    %% steady state computation
    steady_state = 0;
	switch model.advance_type
		case 'solveonestep'
            while ~steady_state
                x_old = x;
                y_old = y;
                % update
				% the most generic one                                                                 
				x = model.solveonestep(x,u_steady,dt);                                          
                y = g(x);
                % increments
                incr_x = norm(x-x_old)/norm(x);
                incr_y = norm((y-y_old)./y_norm);
                if incr_x/dt < opt.eps_x && incr_y/dt < opt.eps_y
                    steady_state = 1;
                end
            end
		case 'solveonestep_linear'
            while ~steady_state
                x_old = x;
                y_old = y;
                % update
				% like "solveonestep" but this entails the solution of a linear system                 
				[K,rhs] = model.solveonestep_linear(x,u_steady,dt);                             
				x = K \ rhs;                                                                            
                y = g(x);
                % increments
                incr_x = norm(x-x_old)/norm(x);
                incr_y = norm((y-y_old)./y_norm);
                if incr_x/dt < opt.eps_x && incr_y/dt < opt.eps_y
                    steady_state = 1;
                end
            end
		case 'nonlinear_explicit'
            while ~steady_state
                x_old = x;
                y_old = y;
                % update
				% explicit euler, with generic right-hand side (NB: model.u_implicit should be false!) 
				x = x + dt * model.f(x,u_steady);                                               
                y = g(x);
                % increments
                incr_x = norm(x-x_old)/norm(x);
                incr_y = norm((y-y_old)./y_norm);
                if incr_x/dt < opt.eps_x && incr_y/dt < opt.eps_y
                    steady_state = 1;
                end
            end
		case 'ANN'
            while ~steady_state
                x_old = x;
                y_old = y;
                % update
				% explicit euler, with ANN right-hand side (NB: model.u_implicit should be false!)     
				if nU > 0                                                                                
					x = x + dt * model.ANN([u_steady;x]);                                      
				else                                                                                     
					x = x + dt * model.ANN(x);                                                     
				end                                                                                      
                y = g(x);
                % increments
                incr_x = norm(x-x_old)/norm(x);
                incr_y = norm((y-y_old)./y_norm);
                if incr_x/dt < opt.eps_x && incr_y/dt < opt.eps_y
                    steady_state = 1;
                end
            end
		case 'linear_CN_uaffine'
            while ~steady_state
                x_old = x;
                y_old = y;
                % update
			   % generic linear model, with Crank-Nicolson, with affine dependence on u                
			   % dx/dt = A x + (b0 + sum_j b_j u_j)                                                    
			   u = u_steady;                 
			   b  = model.b0;                                                                           
			   for iU = 1:nU                                                                            
				   b  = b  + u(iU)*model.b{iU};                                                         
			   end                                                                                      
			   x = (eye(nX)/dt - .5*model.A) \ ((eye(nX)/dt + .5*model.A)*x + b);            
                y = g(x);
                % increments
                incr_x = norm(x-x_old)/norm(x);
                incr_y = norm((y-y_old)./y_norm);
                if incr_x/dt < opt.eps_x && incr_y/dt < opt.eps_y
                    steady_state = 1;
                end
            end
		case 'linear_advance'
            while ~steady_state
                x_old = x;
                y_old = y;
                % update
				% generic linear model                                                                 
				[K,B,b] = model.linear_advance(u_steady,dt);                                    
				x = K \ (B*x + b);                                                                      
                y = g(x);
                % increments
                incr_x = norm(x-x_old)/norm(x);
                incr_y = norm((y-y_old)./y_norm);
                if incr_x/dt < opt.eps_x && incr_y/dt < opt.eps_y
                    steady_state = 1;
                end
            end
		case 'linear_advance_timeexplicit'
            while ~steady_state
                x_old = x;
                y_old = y;
                % update
				% generic linear model, with explicit dependence on dt                                 
				[Kf,Kv,Bf,Bv,b] = model.linear_advance_timeexplicit(u_steady);                       
				x = (Kf + Kv/dt) \ ((Bf + Bv/dt)*x + b);                                      
                y = g(x);
                % increments
                incr_x = norm(x-x_old)/norm(x);
                incr_y = norm((y-y_old)./y_norm);
                if incr_x/dt < opt.eps_x && incr_y/dt < opt.eps_y
                    steady_state = 1;
                end
            end
		case 'linear_advance_uaffine'
            while ~steady_state
                x_old = x;
                y_old = y;
                % update
				% generic linear model, with affine dependence on u                                     
				% NB: we suppose that the matrices are computed for a given value of dt!                                                        
				u = u_steady;                                                                        
				K = model.K0;                                                                            
				B = model.B0;                                                                            
				b  = model.b0;                                                                           
				for iU = 1:nU                                                                            
				   K = K + u(iU)*model.K(:,:,iU);                                                       
				   B = B + u(iU)*model.B(:,:,iU);                                                       
				   b  = b  + u(iU)*model.b(:,iU);                                                       
				end                                                                                      
				x = K \ (B*x + b);                                                                      
                y = g(x);
                % increments
                incr_x = norm(x-x_old)/norm(x);
                incr_y = norm((y-y_old)./y_norm);
                if incr_x/dt < opt.eps_x && incr_y/dt < opt.eps_y
                    steady_state = 1;
                end
            end
		case 'linear_advance_timeexplicit_uaffine'
            while ~steady_state
                x_old = x;
                y_old = y;
                % update
				% generic linear model, with explicit dependence on dt, with affine dependence on u     
				u = u_steady;                                                                        
				Kf = model.Kf0;                                                                          
				Kv = model.Kv0;                                                                          
				Bf = model.Bf0;                                                                          
				Bv = model.Bv0;                                                                          
				b  = model.b0;                                                                           
				for iU = 1:nU                                                                            
				%             Kf = Kf + u(iU)*model.Kf(:,:,iU);                                         
				%             Kv = Kv + u(iU)*model.Kv(:,:,iU);                                         
				%             Bf = Bf + u(iU)*model.Bf(:,:,iU);                                         
				%             Bv = Bv + u(iU)*model.Bv(:,:,iU);                                         
				%             b  = b  + u(iU)*model.b(:,iU);                                            
				   Kf = Kf + u(iU)*model.Kf{iU};                                                        
				   Kv = Kv + u(iU)*model.Kv{iU};                                                        
				   Bf = Bf + u(iU)*model.Bf{iU};                                                        
				   Bv = Bv + u(iU)*model.Bv{iU};                                                        
				   b  = b  + u(iU)*model.b{iU};                                                         
				end                                                                                      
				x = (Kf + Kv/dt) \ ((Bf + Bv/dt)*x + b); 
                y = g(x);
                % increments
                incr_x = norm(x-x_old)/norm(x);
                incr_y = norm((y-y_old)./y_norm);
                if incr_x/dt < opt.eps_x && incr_y/dt < opt.eps_y
                    steady_state = 1;
                end
            end
		otherwise
			error('unknown advance type %s',model.advance_type)
	end
end
