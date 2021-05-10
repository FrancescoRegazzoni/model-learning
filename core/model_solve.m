function output = model_solve(test,model,opt)
    % MODEL_SOLVE solves a given test with a given model.
    %    test.tt: time (vector, or just extrema)
    %    test.uu: input (vector or handle function)
    %    test.yy: exact output (vector, optional)
    %
    % NB: this function is generated automatically by the script
    % "build_model_solve.m", for efficiency reasons. Do not modify it
    % directly. To modify it:
    %  - enter into folder model_solve_builder;
    %  - modifiy the desired section (model_solve_header, model_solve_footer,
    %    etc.);
    %  - run "build_model_solve.m"
    %  - copy the script "model_solve.m" in this folder into "core"
    
    
    opt.dummy = 0;
    
    if model.problem.metaproblem && ~model.problem.particularized
        error('The associated problem is a non-particularized metaproblem!')
    end
    
    if length(test) > 1
        for iS = 1:length(test)
            output{iS} = model_solve(test{iS},model,opt);
        end
        return
    end
    
    if ~isfield(model,'blackbox')
        model.blackbox = 0;
    end    
    
    %% setting default parameters
    if ~isfield(opt,'save_x')
        opt.save_x = 0;
    end
    if ~isfield(opt,'error_compute')
        opt.error_compute = 0;
    end
    if ~isfield(opt,'save_x_freq')
        opt.save_x_freq = 1;
    end
    if ~isfield(opt,'save_y_ex')
        opt.save_y_ex = 1;
    end
    if ~isfield(opt,'verbose')
        opt.verbose = 0;
    end
    if ~isfield(opt,'do_plot')
        opt.do_plot = 0;
    end
    if ~isfield(opt,'do_plot_x')
        opt.do_plot_x = 0;
    end
    if ~isfield(opt,'interpolation_mode_u')
        opt.interpolation_mode_u = 'pointwise';
%         opt.interpolation_mode_u = 'mean_forward';
    end
    if ~isfield(opt,'interpolation_mode_y')
        opt.interpolation_mode_y = 'pointwise';
    end
    
    
    if model.blackbox
        
        init = tic();          
        output = model.forward_function(test,opt);
        timeElapsed = toc(init);
        
        if ~isfield(output,'time')            
            output.time = timeElapsed;
        end
        if ~isfield(output,'time_norm')            
            output.time_norm = timeElapsed/(output.tt(end)-output.tt(1));
        end
        if ~isfield(output,'tt_y')            
            output.tt_y = output.tt;
        end
        
    else

        %% custom initialization
        switch model.advance_type
            case 'nonlinear_explicit'
                if isfield(model,'u_implicit')
                    if model.u_implicit
                        error('advance_type requires u_implicit = 0')
                    end
                end
                model.u_implicit = 0;
        end
        
        %% building time vector
        tt = test.tt(1):model.dt:test.tt(end);
        dtt = [0 tt(2:end)-tt(1:end-1)];

        if ~isfield(test,'tt_y')
            test.tt_y = test.tt;
        end
        idx_tt_y_first = find(tt >= test.tt_y(1),1,'first');
        idx_tt_y_last  = find(tt <= test.tt_y(end),1,'last');
        tt_y = tt(idx_tt_y_first:idx_tt_y_last);

        %% setting constants
        nT = length(tt);
        nT_y = length(tt_y);
        nX = model.nX;
        nY = model.problem.nY;
        nU = model.problem.nU;

        if strcmp(model.advance_type,'ANN')
            if model.nX == nY
                opt.do_plot_x = 0;
            end
        end

        %% building input/output vectors
        if nU > 0        
            if isa(test.uu,'function_handle')
                ufunc = test.uu;
                uu = ufunc(tt);
            else
                uu = interp_time_series(test.tt,test.uu,tt,struct('mode',opt.interpolation_mode_u));
%                 if isequal(tt,test.tt)
%                     uu = test.uu;
%                 else
%                     for iU = 1:size(test.uu,1)
%                         uu(iU,:) = interp1(test.tt,test.uu(iU,:),tt);
%                     end
%                 end
            end

            if model.u_implicit
                uu_eff = uu;
            else
                uu_eff = [zeros(nU,1) uu(:,1:end-1)];
            end
        else
            % I set it to zero, otherwise I get errors while extracting the current u
            uu_eff = zeros(1,nT);
        end

        if isfield(test,'yy') && (opt.error_compute || opt.do_plot || opt.save_y_ex)
            yy_ex = interp_time_series(test.tt_y,test.yy,tt_y,struct('mode',opt.interpolation_mode_y));
%             if isequal(tt_y,test.tt_y)   
%                 yy_ex = test.yy; 
%             else         
%                 for iY = 1:size(test.yy,1)
%                     yy_ex(iY,:) = interp1(test.tt_y,test.yy(iY,:),tt_y);
%                 end
%             end        
        end

        %% initialiation 
        x = model.x0;

        if opt.save_x || opt.do_plot_x
            xx = x;
        end

        yy = model_get_output(model,x);
        
        %% time loop  
        timeInit = tic();
		switch model.output_type
		case 'insidestate'
			if opt.do_plot_x || (opt.save_x && opt.save_x_freq == 1 )
				switch model.advance_type
					case 'solveonestep'
						for iT = 2:nT
							% the most generic one                                                                 
							x = model.solveonestep(x,uu_eff(:,iT),dtt(iT));                                          
							if iT >= idx_tt_y_first && iT<= idx_tt_y_last
								yy(:,iT - idx_tt_y_first + 1) = x(1:nY);
							end
							xx  = [xx x];
						end
					case 'solveonestep_linear'
						for iT = 2:nT
							% like "solveonestep" but this entails the solution of a linear system                 
							[K,rhs] = model.solveonestep_linear(x,uu_eff(:,iT),dtt(iT));                             
							x = K \ rhs;                                                                            
							if iT >= idx_tt_y_first && iT<= idx_tt_y_last
								yy(:,iT - idx_tt_y_first + 1) = x(1:nY);
							end
							xx  = [xx x];
						end
					case 'nonlinear_explicit'
						for iT = 2:nT
							% explicit euler, with generic right-hand side (NB: model.u_implicit should be false!) 
							x = x + dtt(iT) * model.f(x,uu_eff(:,iT));                                               
							if iT >= idx_tt_y_first && iT<= idx_tt_y_last
								yy(:,iT - idx_tt_y_first + 1) = x(1:nY);
							end
							xx  = [xx x];
						end
					case 'ANN'
						for iT = 2:nT
							% explicit euler, with ANN right-hand side (NB: model.u_implicit should be false!)     
							if nU > 0                                                                                
								x = x + dtt(iT) * model.ANN([uu_eff(:,iT);x]);                                      
							else                                                                                     
								x = x + dtt(iT) * model.ANN(x);                                                     
							end                                                                                      
							if iT >= idx_tt_y_first && iT<= idx_tt_y_last
								yy(:,iT - idx_tt_y_first + 1) = x(1:nY);
							end
							xx  = [xx x];
						end
					case 'linear_CN_uaffine'
						for iT = 2:nT
						   % generic linear model, with Crank-Nicolson, with affine dependence on u                
						   % dx/dt = A x + (b0 + sum_j b_j u_j)                                                    
						   u = uu_eff(:,iT);                 
						   b  = model.b0;                                                                           
						   for iU = 1:nU                                                                            
							   b  = b  + u(iU)*model.b{iU};                                                         
						   end                                                                                      
						   x = (eye(nX)/dtt(iT) - .5*model.A) \ ((eye(nX)/dtt(iT) + .5*model.A)*x + b);            
							if iT >= idx_tt_y_first && iT<= idx_tt_y_last
								yy(:,iT - idx_tt_y_first + 1) = x(1:nY);
							end
							xx  = [xx x];
						end
					case 'linear_advance'
						for iT = 2:nT
							% generic linear model                                                                 
							[K,B,b] = model.linear_advance(uu_eff(:,iT),dtt(iT));                                    
							x = K \ (B*x + b);                                                                      
							if iT >= idx_tt_y_first && iT<= idx_tt_y_last
								yy(:,iT - idx_tt_y_first + 1) = x(1:nY);
							end
							xx  = [xx x];
						end
					case 'linear_advance_timeexplicit'
						for iT = 2:nT
							% generic linear model, with explicit dependence on dt                                 
							[Kf,Kv,Bf,Bv,b] = model.linear_advance_timeexplicit(uu_eff(:,iT));                       
							x = (Kf + Kv/dtt(iT)) \ ((Bf + Bv/dtt(iT))*x + b);                                      
							if iT >= idx_tt_y_first && iT<= idx_tt_y_last
								yy(:,iT - idx_tt_y_first + 1) = x(1:nY);
							end
							xx  = [xx x];
						end
					case 'linear_advance_uaffine'
						for iT = 2:nT
							% generic linear model, with affine dependence on u                                     
							% NB: we suppose that the matrices are computed for a given value of dt!                                                        
							u = uu_eff(:,iT);                                                                        
							K = model.K0;                                                                            
							B = model.B0;                                                                            
							b  = model.b0;                                                                           
							for iU = 1:nU                                                                            
							   K = K + u(iU)*model.K(:,:,iU);                                                       
							   B = B + u(iU)*model.B(:,:,iU);                                                       
							   b  = b  + u(iU)*model.b(:,iU);                                                       
							end                                                                                      
							x = K \ (B*x + b);                                                                      
							if iT >= idx_tt_y_first && iT<= idx_tt_y_last
								yy(:,iT - idx_tt_y_first + 1) = x(1:nY);
							end
							xx  = [xx x];
						end
					case 'linear_advance_timeexplicit_uaffine'
						for iT = 2:nT
							% generic linear model, with explicit dependence on dt, with affine dependence on u     
							u = uu_eff(:,iT);                                                                        
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
							x = (Kf + Kv/dtt(iT)) \ ((Bf + Bv/dtt(iT))*x + b); 
							if iT >= idx_tt_y_first && iT<= idx_tt_y_last
								yy(:,iT - idx_tt_y_first + 1) = x(1:nY);
							end
							xx  = [xx x];
						end
					otherwise
						error('unknown advance type %s',model.advance_type)
				end
			elseif opt.save_x
				switch model.advance_type
					case 'solveonestep'
						for iT = 2:nT
							% the most generic one                                                                 
							x = model.solveonestep(x,uu_eff(:,iT),dtt(iT));                                          
							if iT >= idx_tt_y_first && iT<= idx_tt_y_last
								yy(:,iT - idx_tt_y_first + 1) = x(1:nY);
							end
							if mod(iT,opt.save_x_freq) == 0
							    xx  = [xx x];
							end
						end
					case 'solveonestep_linear'
						for iT = 2:nT
							% like "solveonestep" but this entails the solution of a linear system                 
							[K,rhs] = model.solveonestep_linear(x,uu_eff(:,iT),dtt(iT));                             
							x = K \ rhs;                                                                            
							if iT >= idx_tt_y_first && iT<= idx_tt_y_last
								yy(:,iT - idx_tt_y_first + 1) = x(1:nY);
							end
							if mod(iT,opt.save_x_freq) == 0
							    xx  = [xx x];
							end
						end
					case 'nonlinear_explicit'
						for iT = 2:nT
							% explicit euler, with generic right-hand side (NB: model.u_implicit should be false!) 
							x = x + dtt(iT) * model.f(x,uu_eff(:,iT));                                               
							if iT >= idx_tt_y_first && iT<= idx_tt_y_last
								yy(:,iT - idx_tt_y_first + 1) = x(1:nY);
							end
							if mod(iT,opt.save_x_freq) == 0
							    xx  = [xx x];
							end
						end
					case 'ANN'
						for iT = 2:nT
							% explicit euler, with ANN right-hand side (NB: model.u_implicit should be false!)     
							if nU > 0                                                                                
								x = x + dtt(iT) * model.ANN([uu_eff(:,iT);x]);                                      
							else                                                                                     
								x = x + dtt(iT) * model.ANN(x);                                                     
							end                                                                                      
							if iT >= idx_tt_y_first && iT<= idx_tt_y_last
								yy(:,iT - idx_tt_y_first + 1) = x(1:nY);
							end
							if mod(iT,opt.save_x_freq) == 0
							    xx  = [xx x];
							end
						end
					case 'linear_CN_uaffine'
						for iT = 2:nT
						   % generic linear model, with Crank-Nicolson, with affine dependence on u                
						   % dx/dt = A x + (b0 + sum_j b_j u_j)                                                    
						   u = uu_eff(:,iT);                 
						   b  = model.b0;                                                                           
						   for iU = 1:nU                                                                            
							   b  = b  + u(iU)*model.b{iU};                                                         
						   end                                                                                      
						   x = (eye(nX)/dtt(iT) - .5*model.A) \ ((eye(nX)/dtt(iT) + .5*model.A)*x + b);            
							if iT >= idx_tt_y_first && iT<= idx_tt_y_last
								yy(:,iT - idx_tt_y_first + 1) = x(1:nY);
							end
							if mod(iT,opt.save_x_freq) == 0
							    xx  = [xx x];
							end
						end
					case 'linear_advance'
						for iT = 2:nT
							% generic linear model                                                                 
							[K,B,b] = model.linear_advance(uu_eff(:,iT),dtt(iT));                                    
							x = K \ (B*x + b);                                                                      
							if iT >= idx_tt_y_first && iT<= idx_tt_y_last
								yy(:,iT - idx_tt_y_first + 1) = x(1:nY);
							end
							if mod(iT,opt.save_x_freq) == 0
							    xx  = [xx x];
							end
						end
					case 'linear_advance_timeexplicit'
						for iT = 2:nT
							% generic linear model, with explicit dependence on dt                                 
							[Kf,Kv,Bf,Bv,b] = model.linear_advance_timeexplicit(uu_eff(:,iT));                       
							x = (Kf + Kv/dtt(iT)) \ ((Bf + Bv/dtt(iT))*x + b);                                      
							if iT >= idx_tt_y_first && iT<= idx_tt_y_last
								yy(:,iT - idx_tt_y_first + 1) = x(1:nY);
							end
							if mod(iT,opt.save_x_freq) == 0
							    xx  = [xx x];
							end
						end
					case 'linear_advance_uaffine'
						for iT = 2:nT
							% generic linear model, with affine dependence on u                                     
							% NB: we suppose that the matrices are computed for a given value of dt!                                                        
							u = uu_eff(:,iT);                                                                        
							K = model.K0;                                                                            
							B = model.B0;                                                                            
							b  = model.b0;                                                                           
							for iU = 1:nU                                                                            
							   K = K + u(iU)*model.K(:,:,iU);                                                       
							   B = B + u(iU)*model.B(:,:,iU);                                                       
							   b  = b  + u(iU)*model.b(:,iU);                                                       
							end                                                                                      
							x = K \ (B*x + b);                                                                      
							if iT >= idx_tt_y_first && iT<= idx_tt_y_last
								yy(:,iT - idx_tt_y_first + 1) = x(1:nY);
							end
							if mod(iT,opt.save_x_freq) == 0
							    xx  = [xx x];
							end
						end
					case 'linear_advance_timeexplicit_uaffine'
						for iT = 2:nT
							% generic linear model, with explicit dependence on dt, with affine dependence on u     
							u = uu_eff(:,iT);                                                                        
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
							x = (Kf + Kv/dtt(iT)) \ ((Bf + Bv/dtt(iT))*x + b); 
							if iT >= idx_tt_y_first && iT<= idx_tt_y_last
								yy(:,iT - idx_tt_y_first + 1) = x(1:nY);
							end
							if mod(iT,opt.save_x_freq) == 0
							    xx  = [xx x];
							end
						end
					otherwise
						error('unknown advance type %s',model.advance_type)
				end
			else
				switch model.advance_type
					case 'solveonestep'
						for iT = 2:nT
							% the most generic one                                                                 
							x = model.solveonestep(x,uu_eff(:,iT),dtt(iT));                                          
							if iT >= idx_tt_y_first && iT<= idx_tt_y_last
								yy(:,iT - idx_tt_y_first + 1) = x(1:nY);
							end
						end
					case 'solveonestep_linear'
						for iT = 2:nT
							% like "solveonestep" but this entails the solution of a linear system                 
							[K,rhs] = model.solveonestep_linear(x,uu_eff(:,iT),dtt(iT));                             
							x = K \ rhs;                                                                            
							if iT >= idx_tt_y_first && iT<= idx_tt_y_last
								yy(:,iT - idx_tt_y_first + 1) = x(1:nY);
							end
						end
					case 'nonlinear_explicit'
						for iT = 2:nT
							% explicit euler, with generic right-hand side (NB: model.u_implicit should be false!) 
							x = x + dtt(iT) * model.f(x,uu_eff(:,iT));                                               
							if iT >= idx_tt_y_first && iT<= idx_tt_y_last
								yy(:,iT - idx_tt_y_first + 1) = x(1:nY);
							end
						end
					case 'ANN'
						for iT = 2:nT
							% explicit euler, with ANN right-hand side (NB: model.u_implicit should be false!)     
							if nU > 0                                                                                
								x = x + dtt(iT) * model.ANN([uu_eff(:,iT);x]);                                      
							else                                                                                     
								x = x + dtt(iT) * model.ANN(x);                                                     
							end                                                                                      
							if iT >= idx_tt_y_first && iT<= idx_tt_y_last
								yy(:,iT - idx_tt_y_first + 1) = x(1:nY);
							end
						end
					case 'linear_CN_uaffine'
						for iT = 2:nT
						   % generic linear model, with Crank-Nicolson, with affine dependence on u                
						   % dx/dt = A x + (b0 + sum_j b_j u_j)                                                    
						   u = uu_eff(:,iT);                 
						   b  = model.b0;                                                                           
						   for iU = 1:nU                                                                            
							   b  = b  + u(iU)*model.b{iU};                                                         
						   end                                                                                      
						   x = (eye(nX)/dtt(iT) - .5*model.A) \ ((eye(nX)/dtt(iT) + .5*model.A)*x + b);            
							if iT >= idx_tt_y_first && iT<= idx_tt_y_last
								yy(:,iT - idx_tt_y_first + 1) = x(1:nY);
							end
						end
					case 'linear_advance'
						for iT = 2:nT
							% generic linear model                                                                 
							[K,B,b] = model.linear_advance(uu_eff(:,iT),dtt(iT));                                    
							x = K \ (B*x + b);                                                                      
							if iT >= idx_tt_y_first && iT<= idx_tt_y_last
								yy(:,iT - idx_tt_y_first + 1) = x(1:nY);
							end
						end
					case 'linear_advance_timeexplicit'
						for iT = 2:nT
							% generic linear model, with explicit dependence on dt                                 
							[Kf,Kv,Bf,Bv,b] = model.linear_advance_timeexplicit(uu_eff(:,iT));                       
							x = (Kf + Kv/dtt(iT)) \ ((Bf + Bv/dtt(iT))*x + b);                                      
							if iT >= idx_tt_y_first && iT<= idx_tt_y_last
								yy(:,iT - idx_tt_y_first + 1) = x(1:nY);
							end
						end
					case 'linear_advance_uaffine'
						for iT = 2:nT
							% generic linear model, with affine dependence on u                                     
							% NB: we suppose that the matrices are computed for a given value of dt!                                                        
							u = uu_eff(:,iT);                                                                        
							K = model.K0;                                                                            
							B = model.B0;                                                                            
							b  = model.b0;                                                                           
							for iU = 1:nU                                                                            
							   K = K + u(iU)*model.K(:,:,iU);                                                       
							   B = B + u(iU)*model.B(:,:,iU);                                                       
							   b  = b  + u(iU)*model.b(:,iU);                                                       
							end                                                                                      
							x = K \ (B*x + b);                                                                      
							if iT >= idx_tt_y_first && iT<= idx_tt_y_last
								yy(:,iT - idx_tt_y_first + 1) = x(1:nY);
							end
						end
					case 'linear_advance_timeexplicit_uaffine'
						for iT = 2:nT
							% generic linear model, with explicit dependence on dt, with affine dependence on u     
							u = uu_eff(:,iT);                                                                        
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
							x = (Kf + Kv/dtt(iT)) \ ((Bf + Bv/dtt(iT))*x + b); 
							if iT >= idx_tt_y_first && iT<= idx_tt_y_last
								yy(:,iT - idx_tt_y_first + 1) = x(1:nY);
							end
						end
					otherwise
						error('unknown advance type %s',model.advance_type)
				end
			end
		case 'nonlinear'
			if opt.do_plot_x || (opt.save_x && opt.save_x_freq == 1 )
				switch model.advance_type
					case 'solveonestep'
						for iT = 2:nT
							% the most generic one                                                                 
							x = model.solveonestep(x,uu_eff(:,iT),dtt(iT));                                          
							if iT >= idx_tt_y_first && iT<= idx_tt_y_last
								yy(:,iT - idx_tt_y_first + 1) = model.g(x);
							end
							xx  = [xx x];
						end
					case 'solveonestep_linear'
						for iT = 2:nT
							% like "solveonestep" but this entails the solution of a linear system                 
							[K,rhs] = model.solveonestep_linear(x,uu_eff(:,iT),dtt(iT));                             
							x = K \ rhs;                                                                            
							if iT >= idx_tt_y_first && iT<= idx_tt_y_last
								yy(:,iT - idx_tt_y_first + 1) = model.g(x);
							end
							xx  = [xx x];
						end
					case 'nonlinear_explicit'
						for iT = 2:nT
							% explicit euler, with generic right-hand side (NB: model.u_implicit should be false!) 
							x = x + dtt(iT) * model.f(x,uu_eff(:,iT));                                               
							if iT >= idx_tt_y_first && iT<= idx_tt_y_last
								yy(:,iT - idx_tt_y_first + 1) = model.g(x);
							end
							xx  = [xx x];
						end
					case 'ANN'
						for iT = 2:nT
							% explicit euler, with ANN right-hand side (NB: model.u_implicit should be false!)     
							if nU > 0                                                                                
								x = x + dtt(iT) * model.ANN([uu_eff(:,iT);x]);                                      
							else                                                                                     
								x = x + dtt(iT) * model.ANN(x);                                                     
							end                                                                                      
							if iT >= idx_tt_y_first && iT<= idx_tt_y_last
								yy(:,iT - idx_tt_y_first + 1) = model.g(x);
							end
							xx  = [xx x];
						end
					case 'linear_CN_uaffine'
						for iT = 2:nT
						   % generic linear model, with Crank-Nicolson, with affine dependence on u                
						   % dx/dt = A x + (b0 + sum_j b_j u_j)                                                    
						   u = uu_eff(:,iT);                 
						   b  = model.b0;                                                                           
						   for iU = 1:nU                                                                            
							   b  = b  + u(iU)*model.b{iU};                                                         
						   end                                                                                      
						   x = (eye(nX)/dtt(iT) - .5*model.A) \ ((eye(nX)/dtt(iT) + .5*model.A)*x + b);            
							if iT >= idx_tt_y_first && iT<= idx_tt_y_last
								yy(:,iT - idx_tt_y_first + 1) = model.g(x);
							end
							xx  = [xx x];
						end
					case 'linear_advance'
						for iT = 2:nT
							% generic linear model                                                                 
							[K,B,b] = model.linear_advance(uu_eff(:,iT),dtt(iT));                                    
							x = K \ (B*x + b);                                                                      
							if iT >= idx_tt_y_first && iT<= idx_tt_y_last
								yy(:,iT - idx_tt_y_first + 1) = model.g(x);
							end
							xx  = [xx x];
						end
					case 'linear_advance_timeexplicit'
						for iT = 2:nT
							% generic linear model, with explicit dependence on dt                                 
							[Kf,Kv,Bf,Bv,b] = model.linear_advance_timeexplicit(uu_eff(:,iT));                       
							x = (Kf + Kv/dtt(iT)) \ ((Bf + Bv/dtt(iT))*x + b);                                      
							if iT >= idx_tt_y_first && iT<= idx_tt_y_last
								yy(:,iT - idx_tt_y_first + 1) = model.g(x);
							end
							xx  = [xx x];
						end
					case 'linear_advance_uaffine'
						for iT = 2:nT
							% generic linear model, with affine dependence on u                                     
							% NB: we suppose that the matrices are computed for a given value of dt!                                                        
							u = uu_eff(:,iT);                                                                        
							K = model.K0;                                                                            
							B = model.B0;                                                                            
							b  = model.b0;                                                                           
							for iU = 1:nU                                                                            
							   K = K + u(iU)*model.K(:,:,iU);                                                       
							   B = B + u(iU)*model.B(:,:,iU);                                                       
							   b  = b  + u(iU)*model.b(:,iU);                                                       
							end                                                                                      
							x = K \ (B*x + b);                                                                      
							if iT >= idx_tt_y_first && iT<= idx_tt_y_last
								yy(:,iT - idx_tt_y_first + 1) = model.g(x);
							end
							xx  = [xx x];
						end
					case 'linear_advance_timeexplicit_uaffine'
						for iT = 2:nT
							% generic linear model, with explicit dependence on dt, with affine dependence on u     
							u = uu_eff(:,iT);                                                                        
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
							x = (Kf + Kv/dtt(iT)) \ ((Bf + Bv/dtt(iT))*x + b); 
							if iT >= idx_tt_y_first && iT<= idx_tt_y_last
								yy(:,iT - idx_tt_y_first + 1) = model.g(x);
							end
							xx  = [xx x];
						end
					otherwise
						error('unknown advance type %s',model.advance_type)
				end
			elseif opt.save_x
				switch model.advance_type
					case 'solveonestep'
						for iT = 2:nT
							% the most generic one                                                                 
							x = model.solveonestep(x,uu_eff(:,iT),dtt(iT));                                          
							if iT >= idx_tt_y_first && iT<= idx_tt_y_last
								yy(:,iT - idx_tt_y_first + 1) = model.g(x);
							end
							if mod(iT,opt.save_x_freq) == 0
							    xx  = [xx x];
							end
						end
					case 'solveonestep_linear'
						for iT = 2:nT
							% like "solveonestep" but this entails the solution of a linear system                 
							[K,rhs] = model.solveonestep_linear(x,uu_eff(:,iT),dtt(iT));                             
							x = K \ rhs;                                                                            
							if iT >= idx_tt_y_first && iT<= idx_tt_y_last
								yy(:,iT - idx_tt_y_first + 1) = model.g(x);
							end
							if mod(iT,opt.save_x_freq) == 0
							    xx  = [xx x];
							end
						end
					case 'nonlinear_explicit'
						for iT = 2:nT
							% explicit euler, with generic right-hand side (NB: model.u_implicit should be false!) 
							x = x + dtt(iT) * model.f(x,uu_eff(:,iT));                                               
							if iT >= idx_tt_y_first && iT<= idx_tt_y_last
								yy(:,iT - idx_tt_y_first + 1) = model.g(x);
							end
							if mod(iT,opt.save_x_freq) == 0
							    xx  = [xx x];
							end
						end
					case 'ANN'
						for iT = 2:nT
							% explicit euler, with ANN right-hand side (NB: model.u_implicit should be false!)     
							if nU > 0                                                                                
								x = x + dtt(iT) * model.ANN([uu_eff(:,iT);x]);                                      
							else                                                                                     
								x = x + dtt(iT) * model.ANN(x);                                                     
							end                                                                                      
							if iT >= idx_tt_y_first && iT<= idx_tt_y_last
								yy(:,iT - idx_tt_y_first + 1) = model.g(x);
							end
							if mod(iT,opt.save_x_freq) == 0
							    xx  = [xx x];
							end
						end
					case 'linear_CN_uaffine'
						for iT = 2:nT
						   % generic linear model, with Crank-Nicolson, with affine dependence on u                
						   % dx/dt = A x + (b0 + sum_j b_j u_j)                                                    
						   u = uu_eff(:,iT);                 
						   b  = model.b0;                                                                           
						   for iU = 1:nU                                                                            
							   b  = b  + u(iU)*model.b{iU};                                                         
						   end                                                                                      
						   x = (eye(nX)/dtt(iT) - .5*model.A) \ ((eye(nX)/dtt(iT) + .5*model.A)*x + b);            
							if iT >= idx_tt_y_first && iT<= idx_tt_y_last
								yy(:,iT - idx_tt_y_first + 1) = model.g(x);
							end
							if mod(iT,opt.save_x_freq) == 0
							    xx  = [xx x];
							end
						end
					case 'linear_advance'
						for iT = 2:nT
							% generic linear model                                                                 
							[K,B,b] = model.linear_advance(uu_eff(:,iT),dtt(iT));                                    
							x = K \ (B*x + b);                                                                      
							if iT >= idx_tt_y_first && iT<= idx_tt_y_last
								yy(:,iT - idx_tt_y_first + 1) = model.g(x);
							end
							if mod(iT,opt.save_x_freq) == 0
							    xx  = [xx x];
							end
						end
					case 'linear_advance_timeexplicit'
						for iT = 2:nT
							% generic linear model, with explicit dependence on dt                                 
							[Kf,Kv,Bf,Bv,b] = model.linear_advance_timeexplicit(uu_eff(:,iT));                       
							x = (Kf + Kv/dtt(iT)) \ ((Bf + Bv/dtt(iT))*x + b);                                      
							if iT >= idx_tt_y_first && iT<= idx_tt_y_last
								yy(:,iT - idx_tt_y_first + 1) = model.g(x);
							end
							if mod(iT,opt.save_x_freq) == 0
							    xx  = [xx x];
							end
						end
					case 'linear_advance_uaffine'
						for iT = 2:nT
							% generic linear model, with affine dependence on u                                     
							% NB: we suppose that the matrices are computed for a given value of dt!                                                        
							u = uu_eff(:,iT);                                                                        
							K = model.K0;                                                                            
							B = model.B0;                                                                            
							b  = model.b0;                                                                           
							for iU = 1:nU                                                                            
							   K = K + u(iU)*model.K(:,:,iU);                                                       
							   B = B + u(iU)*model.B(:,:,iU);                                                       
							   b  = b  + u(iU)*model.b(:,iU);                                                       
							end                                                                                      
							x = K \ (B*x + b);                                                                      
							if iT >= idx_tt_y_first && iT<= idx_tt_y_last
								yy(:,iT - idx_tt_y_first + 1) = model.g(x);
							end
							if mod(iT,opt.save_x_freq) == 0
							    xx  = [xx x];
							end
						end
					case 'linear_advance_timeexplicit_uaffine'
						for iT = 2:nT
							% generic linear model, with explicit dependence on dt, with affine dependence on u     
							u = uu_eff(:,iT);                                                                        
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
							x = (Kf + Kv/dtt(iT)) \ ((Bf + Bv/dtt(iT))*x + b); 
							if iT >= idx_tt_y_first && iT<= idx_tt_y_last
								yy(:,iT - idx_tt_y_first + 1) = model.g(x);
							end
							if mod(iT,opt.save_x_freq) == 0
							    xx  = [xx x];
							end
						end
					otherwise
						error('unknown advance type %s',model.advance_type)
				end
			else
				switch model.advance_type
					case 'solveonestep'
						for iT = 2:nT
							% the most generic one                                                                 
							x = model.solveonestep(x,uu_eff(:,iT),dtt(iT));                                          
							if iT >= idx_tt_y_first && iT<= idx_tt_y_last
								yy(:,iT - idx_tt_y_first + 1) = model.g(x);
							end
						end
					case 'solveonestep_linear'
						for iT = 2:nT
							% like "solveonestep" but this entails the solution of a linear system                 
							[K,rhs] = model.solveonestep_linear(x,uu_eff(:,iT),dtt(iT));                             
							x = K \ rhs;                                                                            
							if iT >= idx_tt_y_first && iT<= idx_tt_y_last
								yy(:,iT - idx_tt_y_first + 1) = model.g(x);
							end
						end
					case 'nonlinear_explicit'
						for iT = 2:nT
							% explicit euler, with generic right-hand side (NB: model.u_implicit should be false!) 
							x = x + dtt(iT) * model.f(x,uu_eff(:,iT));                                               
							if iT >= idx_tt_y_first && iT<= idx_tt_y_last
								yy(:,iT - idx_tt_y_first + 1) = model.g(x);
							end
						end
					case 'ANN'
						for iT = 2:nT
							% explicit euler, with ANN right-hand side (NB: model.u_implicit should be false!)     
							if nU > 0                                                                                
								x = x + dtt(iT) * model.ANN([uu_eff(:,iT);x]);                                      
							else                                                                                     
								x = x + dtt(iT) * model.ANN(x);                                                     
							end                                                                                      
							if iT >= idx_tt_y_first && iT<= idx_tt_y_last
								yy(:,iT - idx_tt_y_first + 1) = model.g(x);
							end
						end
					case 'linear_CN_uaffine'
						for iT = 2:nT
						   % generic linear model, with Crank-Nicolson, with affine dependence on u                
						   % dx/dt = A x + (b0 + sum_j b_j u_j)                                                    
						   u = uu_eff(:,iT);                 
						   b  = model.b0;                                                                           
						   for iU = 1:nU                                                                            
							   b  = b  + u(iU)*model.b{iU};                                                         
						   end                                                                                      
						   x = (eye(nX)/dtt(iT) - .5*model.A) \ ((eye(nX)/dtt(iT) + .5*model.A)*x + b);            
							if iT >= idx_tt_y_first && iT<= idx_tt_y_last
								yy(:,iT - idx_tt_y_first + 1) = model.g(x);
							end
						end
					case 'linear_advance'
						for iT = 2:nT
							% generic linear model                                                                 
							[K,B,b] = model.linear_advance(uu_eff(:,iT),dtt(iT));                                    
							x = K \ (B*x + b);                                                                      
							if iT >= idx_tt_y_first && iT<= idx_tt_y_last
								yy(:,iT - idx_tt_y_first + 1) = model.g(x);
							end
						end
					case 'linear_advance_timeexplicit'
						for iT = 2:nT
							% generic linear model, with explicit dependence on dt                                 
							[Kf,Kv,Bf,Bv,b] = model.linear_advance_timeexplicit(uu_eff(:,iT));                       
							x = (Kf + Kv/dtt(iT)) \ ((Bf + Bv/dtt(iT))*x + b);                                      
							if iT >= idx_tt_y_first && iT<= idx_tt_y_last
								yy(:,iT - idx_tt_y_first + 1) = model.g(x);
							end
						end
					case 'linear_advance_uaffine'
						for iT = 2:nT
							% generic linear model, with affine dependence on u                                     
							% NB: we suppose that the matrices are computed for a given value of dt!                                                        
							u = uu_eff(:,iT);                                                                        
							K = model.K0;                                                                            
							B = model.B0;                                                                            
							b  = model.b0;                                                                           
							for iU = 1:nU                                                                            
							   K = K + u(iU)*model.K(:,:,iU);                                                       
							   B = B + u(iU)*model.B(:,:,iU);                                                       
							   b  = b  + u(iU)*model.b(:,iU);                                                       
							end                                                                                      
							x = K \ (B*x + b);                                                                      
							if iT >= idx_tt_y_first && iT<= idx_tt_y_last
								yy(:,iT - idx_tt_y_first + 1) = model.g(x);
							end
						end
					case 'linear_advance_timeexplicit_uaffine'
						for iT = 2:nT
							% generic linear model, with explicit dependence on dt, with affine dependence on u     
							u = uu_eff(:,iT);                                                                        
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
							x = (Kf + Kv/dtt(iT)) \ ((Bf + Bv/dtt(iT))*x + b); 
							if iT >= idx_tt_y_first && iT<= idx_tt_y_last
								yy(:,iT - idx_tt_y_first + 1) = model.g(x);
							end
						end
					otherwise
						error('unknown advance type %s',model.advance_type)
				end
			end
		case 'linear'
			if opt.do_plot_x || (opt.save_x && opt.save_x_freq == 1 )
				switch model.advance_type
					case 'solveonestep'
						for iT = 2:nT
							% the most generic one                                                                 
							x = model.solveonestep(x,uu_eff(:,iT),dtt(iT));                                          
							if iT >= idx_tt_y_first && iT<= idx_tt_y_last
								yy(:,iT - idx_tt_y_first + 1) = model.G*x+model.g0;
							end
							xx  = [xx x];
						end
					case 'solveonestep_linear'
						for iT = 2:nT
							% like "solveonestep" but this entails the solution of a linear system                 
							[K,rhs] = model.solveonestep_linear(x,uu_eff(:,iT),dtt(iT));                             
							x = K \ rhs;                                                                            
							if iT >= idx_tt_y_first && iT<= idx_tt_y_last
								yy(:,iT - idx_tt_y_first + 1) = model.G*x+model.g0;
							end
							xx  = [xx x];
						end
					case 'nonlinear_explicit'
						for iT = 2:nT
							% explicit euler, with generic right-hand side (NB: model.u_implicit should be false!) 
							x = x + dtt(iT) * model.f(x,uu_eff(:,iT));                                               
							if iT >= idx_tt_y_first && iT<= idx_tt_y_last
								yy(:,iT - idx_tt_y_first + 1) = model.G*x+model.g0;
							end
							xx  = [xx x];
						end
					case 'ANN'
						for iT = 2:nT
							% explicit euler, with ANN right-hand side (NB: model.u_implicit should be false!)     
							if nU > 0                                                                                
								x = x + dtt(iT) * model.ANN([uu_eff(:,iT);x]);                                      
							else                                                                                     
								x = x + dtt(iT) * model.ANN(x);                                                     
							end                                                                                      
							if iT >= idx_tt_y_first && iT<= idx_tt_y_last
								yy(:,iT - idx_tt_y_first + 1) = model.G*x+model.g0;
							end
							xx  = [xx x];
						end
					case 'linear_CN_uaffine'
						for iT = 2:nT
						   % generic linear model, with Crank-Nicolson, with affine dependence on u                
						   % dx/dt = A x + (b0 + sum_j b_j u_j)                                                    
						   u = uu_eff(:,iT);                 
						   b  = model.b0;                                                                           
						   for iU = 1:nU                                                                            
							   b  = b  + u(iU)*model.b{iU};                                                         
						   end                                                                                      
						   x = (eye(nX)/dtt(iT) - .5*model.A) \ ((eye(nX)/dtt(iT) + .5*model.A)*x + b);            
							if iT >= idx_tt_y_first && iT<= idx_tt_y_last
								yy(:,iT - idx_tt_y_first + 1) = model.G*x+model.g0;
							end
							xx  = [xx x];
						end
					case 'linear_advance'
						for iT = 2:nT
							% generic linear model                                                                 
							[K,B,b] = model.linear_advance(uu_eff(:,iT),dtt(iT));                                    
							x = K \ (B*x + b);                                                                      
							if iT >= idx_tt_y_first && iT<= idx_tt_y_last
								yy(:,iT - idx_tt_y_first + 1) = model.G*x+model.g0;
							end
							xx  = [xx x];
						end
					case 'linear_advance_timeexplicit'
						for iT = 2:nT
							% generic linear model, with explicit dependence on dt                                 
							[Kf,Kv,Bf,Bv,b] = model.linear_advance_timeexplicit(uu_eff(:,iT));                       
							x = (Kf + Kv/dtt(iT)) \ ((Bf + Bv/dtt(iT))*x + b);                                      
							if iT >= idx_tt_y_first && iT<= idx_tt_y_last
								yy(:,iT - idx_tt_y_first + 1) = model.G*x+model.g0;
							end
							xx  = [xx x];
						end
					case 'linear_advance_uaffine'
						for iT = 2:nT
							% generic linear model, with affine dependence on u                                     
							% NB: we suppose that the matrices are computed for a given value of dt!                                                        
							u = uu_eff(:,iT);                                                                        
							K = model.K0;                                                                            
							B = model.B0;                                                                            
							b  = model.b0;                                                                           
							for iU = 1:nU                                                                            
							   K = K + u(iU)*model.K(:,:,iU);                                                       
							   B = B + u(iU)*model.B(:,:,iU);                                                       
							   b  = b  + u(iU)*model.b(:,iU);                                                       
							end                                                                                      
							x = K \ (B*x + b);                                                                      
							if iT >= idx_tt_y_first && iT<= idx_tt_y_last
								yy(:,iT - idx_tt_y_first + 1) = model.G*x+model.g0;
							end
							xx  = [xx x];
						end
					case 'linear_advance_timeexplicit_uaffine'
						for iT = 2:nT
							% generic linear model, with explicit dependence on dt, with affine dependence on u     
							u = uu_eff(:,iT);                                                                        
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
							x = (Kf + Kv/dtt(iT)) \ ((Bf + Bv/dtt(iT))*x + b); 
							if iT >= idx_tt_y_first && iT<= idx_tt_y_last
								yy(:,iT - idx_tt_y_first + 1) = model.G*x+model.g0;
							end
							xx  = [xx x];
						end
					otherwise
						error('unknown advance type %s',model.advance_type)
				end
			elseif opt.save_x
				switch model.advance_type
					case 'solveonestep'
						for iT = 2:nT
							% the most generic one                                                                 
							x = model.solveonestep(x,uu_eff(:,iT),dtt(iT));                                          
							if iT >= idx_tt_y_first && iT<= idx_tt_y_last
								yy(:,iT - idx_tt_y_first + 1) = model.G*x+model.g0;
							end
							if mod(iT,opt.save_x_freq) == 0
							    xx  = [xx x];
							end
						end
					case 'solveonestep_linear'
						for iT = 2:nT
							% like "solveonestep" but this entails the solution of a linear system                 
							[K,rhs] = model.solveonestep_linear(x,uu_eff(:,iT),dtt(iT));                             
							x = K \ rhs;                                                                            
							if iT >= idx_tt_y_first && iT<= idx_tt_y_last
								yy(:,iT - idx_tt_y_first + 1) = model.G*x+model.g0;
							end
							if mod(iT,opt.save_x_freq) == 0
							    xx  = [xx x];
							end
						end
					case 'nonlinear_explicit'
						for iT = 2:nT
							% explicit euler, with generic right-hand side (NB: model.u_implicit should be false!) 
							x = x + dtt(iT) * model.f(x,uu_eff(:,iT));                                               
							if iT >= idx_tt_y_first && iT<= idx_tt_y_last
								yy(:,iT - idx_tt_y_first + 1) = model.G*x+model.g0;
							end
							if mod(iT,opt.save_x_freq) == 0
							    xx  = [xx x];
							end
						end
					case 'ANN'
						for iT = 2:nT
							% explicit euler, with ANN right-hand side (NB: model.u_implicit should be false!)     
							if nU > 0                                                                                
								x = x + dtt(iT) * model.ANN([uu_eff(:,iT);x]);                                      
							else                                                                                     
								x = x + dtt(iT) * model.ANN(x);                                                     
							end                                                                                      
							if iT >= idx_tt_y_first && iT<= idx_tt_y_last
								yy(:,iT - idx_tt_y_first + 1) = model.G*x+model.g0;
							end
							if mod(iT,opt.save_x_freq) == 0
							    xx  = [xx x];
							end
						end
					case 'linear_CN_uaffine'
						for iT = 2:nT
						   % generic linear model, with Crank-Nicolson, with affine dependence on u                
						   % dx/dt = A x + (b0 + sum_j b_j u_j)                                                    
						   u = uu_eff(:,iT);                 
						   b  = model.b0;                                                                           
						   for iU = 1:nU                                                                            
							   b  = b  + u(iU)*model.b{iU};                                                         
						   end                                                                                      
						   x = (eye(nX)/dtt(iT) - .5*model.A) \ ((eye(nX)/dtt(iT) + .5*model.A)*x + b);            
							if iT >= idx_tt_y_first && iT<= idx_tt_y_last
								yy(:,iT - idx_tt_y_first + 1) = model.G*x+model.g0;
							end
							if mod(iT,opt.save_x_freq) == 0
							    xx  = [xx x];
							end
						end
					case 'linear_advance'
						for iT = 2:nT
							% generic linear model                                                                 
							[K,B,b] = model.linear_advance(uu_eff(:,iT),dtt(iT));                                    
							x = K \ (B*x + b);                                                                      
							if iT >= idx_tt_y_first && iT<= idx_tt_y_last
								yy(:,iT - idx_tt_y_first + 1) = model.G*x+model.g0;
							end
							if mod(iT,opt.save_x_freq) == 0
							    xx  = [xx x];
							end
						end
					case 'linear_advance_timeexplicit'
						for iT = 2:nT
							% generic linear model, with explicit dependence on dt                                 
							[Kf,Kv,Bf,Bv,b] = model.linear_advance_timeexplicit(uu_eff(:,iT));                       
							x = (Kf + Kv/dtt(iT)) \ ((Bf + Bv/dtt(iT))*x + b);                                      
							if iT >= idx_tt_y_first && iT<= idx_tt_y_last
								yy(:,iT - idx_tt_y_first + 1) = model.G*x+model.g0;
							end
							if mod(iT,opt.save_x_freq) == 0
							    xx  = [xx x];
							end
						end
					case 'linear_advance_uaffine'
						for iT = 2:nT
							% generic linear model, with affine dependence on u                                     
							% NB: we suppose that the matrices are computed for a given value of dt!                                                        
							u = uu_eff(:,iT);                                                                        
							K = model.K0;                                                                            
							B = model.B0;                                                                            
							b  = model.b0;                                                                           
							for iU = 1:nU                                                                            
							   K = K + u(iU)*model.K(:,:,iU);                                                       
							   B = B + u(iU)*model.B(:,:,iU);                                                       
							   b  = b  + u(iU)*model.b(:,iU);                                                       
							end                                                                                      
							x = K \ (B*x + b);                                                                      
							if iT >= idx_tt_y_first && iT<= idx_tt_y_last
								yy(:,iT - idx_tt_y_first + 1) = model.G*x+model.g0;
							end
							if mod(iT,opt.save_x_freq) == 0
							    xx  = [xx x];
							end
						end
					case 'linear_advance_timeexplicit_uaffine'
						for iT = 2:nT
							% generic linear model, with explicit dependence on dt, with affine dependence on u     
							u = uu_eff(:,iT);                                                                        
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
							x = (Kf + Kv/dtt(iT)) \ ((Bf + Bv/dtt(iT))*x + b); 
							if iT >= idx_tt_y_first && iT<= idx_tt_y_last
								yy(:,iT - idx_tt_y_first + 1) = model.G*x+model.g0;
							end
							if mod(iT,opt.save_x_freq) == 0
							    xx  = [xx x];
							end
						end
					otherwise
						error('unknown advance type %s',model.advance_type)
				end
			else
				switch model.advance_type
					case 'solveonestep'
						for iT = 2:nT
							% the most generic one                                                                 
							x = model.solveonestep(x,uu_eff(:,iT),dtt(iT));                                          
							if iT >= idx_tt_y_first && iT<= idx_tt_y_last
								yy(:,iT - idx_tt_y_first + 1) = model.G*x+model.g0;
							end
						end
					case 'solveonestep_linear'
						for iT = 2:nT
							% like "solveonestep" but this entails the solution of a linear system                 
							[K,rhs] = model.solveonestep_linear(x,uu_eff(:,iT),dtt(iT));                             
							x = K \ rhs;                                                                            
							if iT >= idx_tt_y_first && iT<= idx_tt_y_last
								yy(:,iT - idx_tt_y_first + 1) = model.G*x+model.g0;
							end
						end
					case 'nonlinear_explicit'
						for iT = 2:nT
							% explicit euler, with generic right-hand side (NB: model.u_implicit should be false!) 
							x = x + dtt(iT) * model.f(x,uu_eff(:,iT));                                               
							if iT >= idx_tt_y_first && iT<= idx_tt_y_last
								yy(:,iT - idx_tt_y_first + 1) = model.G*x+model.g0;
							end
						end
					case 'ANN'
						for iT = 2:nT
							% explicit euler, with ANN right-hand side (NB: model.u_implicit should be false!)     
							if nU > 0                                                                                
								x = x + dtt(iT) * model.ANN([uu_eff(:,iT);x]);                                      
							else                                                                                     
								x = x + dtt(iT) * model.ANN(x);                                                     
							end                                                                                      
							if iT >= idx_tt_y_first && iT<= idx_tt_y_last
								yy(:,iT - idx_tt_y_first + 1) = model.G*x+model.g0;
							end
						end
					case 'linear_CN_uaffine'
						for iT = 2:nT
						   % generic linear model, with Crank-Nicolson, with affine dependence on u                
						   % dx/dt = A x + (b0 + sum_j b_j u_j)                                                    
						   u = uu_eff(:,iT);                 
						   b  = model.b0;                                                                           
						   for iU = 1:nU                                                                            
							   b  = b  + u(iU)*model.b{iU};                                                         
						   end                                                                                      
						   x = (eye(nX)/dtt(iT) - .5*model.A) \ ((eye(nX)/dtt(iT) + .5*model.A)*x + b);            
							if iT >= idx_tt_y_first && iT<= idx_tt_y_last
								yy(:,iT - idx_tt_y_first + 1) = model.G*x+model.g0;
							end
						end
					case 'linear_advance'
						for iT = 2:nT
							% generic linear model                                                                 
							[K,B,b] = model.linear_advance(uu_eff(:,iT),dtt(iT));                                    
							x = K \ (B*x + b);                                                                      
							if iT >= idx_tt_y_first && iT<= idx_tt_y_last
								yy(:,iT - idx_tt_y_first + 1) = model.G*x+model.g0;
							end
						end
					case 'linear_advance_timeexplicit'
						for iT = 2:nT
							% generic linear model, with explicit dependence on dt                                 
							[Kf,Kv,Bf,Bv,b] = model.linear_advance_timeexplicit(uu_eff(:,iT));                       
							x = (Kf + Kv/dtt(iT)) \ ((Bf + Bv/dtt(iT))*x + b);                                      
							if iT >= idx_tt_y_first && iT<= idx_tt_y_last
								yy(:,iT - idx_tt_y_first + 1) = model.G*x+model.g0;
							end
						end
					case 'linear_advance_uaffine'
						for iT = 2:nT
							% generic linear model, with affine dependence on u                                     
							% NB: we suppose that the matrices are computed for a given value of dt!                                                        
							u = uu_eff(:,iT);                                                                        
							K = model.K0;                                                                            
							B = model.B0;                                                                            
							b  = model.b0;                                                                           
							for iU = 1:nU                                                                            
							   K = K + u(iU)*model.K(:,:,iU);                                                       
							   B = B + u(iU)*model.B(:,:,iU);                                                       
							   b  = b  + u(iU)*model.b(:,iU);                                                       
							end                                                                                      
							x = K \ (B*x + b);                                                                      
							if iT >= idx_tt_y_first && iT<= idx_tt_y_last
								yy(:,iT - idx_tt_y_first + 1) = model.G*x+model.g0;
							end
						end
					case 'linear_advance_timeexplicit_uaffine'
						for iT = 2:nT
							% generic linear model, with explicit dependence on dt, with affine dependence on u     
							u = uu_eff(:,iT);                                                                        
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
							x = (Kf + Kv/dtt(iT)) \ ((Bf + Bv/dtt(iT))*x + b); 
							if iT >= idx_tt_y_first && iT<= idx_tt_y_last
								yy(:,iT - idx_tt_y_first + 1) = model.G*x+model.g0;
							end
						end
					otherwise
						error('unknown advance type %s',model.advance_type)
				end
			end
		otherwise
			error('unknown output type %s',model.output_type)
		end
        timeElapsed = toc(timeInit);

        %% postprocessing
        output.tt = tt;
        if nU > 0
            output.uu = uu;
        end
        output.tt_y = tt_y;
        output.yy = yy;
        output.time = timeElapsed;
        output.time_norm = timeElapsed/(tt(end)-tt(1));
        if opt.save_x
            output.xx = xx;
        end
        if opt.error_compute
            output.err_L2 = get_norm_L2_time(tt_y,yy-yy_ex);
    %         output.err_L2 = sqrt(mean((yy(:) - yy_ex(:)).^2));
        end
        if isfield(test,'yy') && opt.save_y_ex
            output.yy_ex = yy_ex;
        end
        if isfield(test,'label')
            output.label = test.label;
        end

        if opt.verbose
            fprintf('   model solved --- time elapsed: %1.2e s', timeElapsed)
            if opt.error_compute
                fprintf(' --- error %1.2e',output.err_L2);
            end
            fprintf('\n');
        end

    end
                    
    if opt.do_plot
        nrows = 1;
        if model.problem.nU > 0
            nrows = nrows+1;
        end
        if opt.do_plot_x 
            nrows = nrows+1;
        end
        
        i_row = 1;
        
        if model.problem.nU > 0
            subplot(nrows,1,i_row)
            plot(output.tt,output.uu,'-','linewidth',1.2)
            axis([0 output.tt(end) min(model.problem.u_min) max(model.problem.u_max)])
            ylabel('u')
            i_row = i_row + 1;
        end
        
        subplot(nrows,1,i_row)
        plot(output.tt_y,output.yy,'-','linewidth',1.2) 
        axis([0 output.tt(end) min(model.problem.y_min) max(model.problem.y_max)])
        if isfield(test,'yy')
            hold on
            sty = '--'; lnwdt = 1.2;
            if nY > 1
                % if the number of outputs ig greater than 1, I use the same colour
                ax = gca;
                ax.ColorOrderIndex = 1;
                sty = '--'; lnwdt = 1.5;
            end
            plot(tt_y,yy_ex,sty,'linewidth',lnwdt)
            if nY == 1
                legend('test model','HF model')
            end
            hold off
        end
        ylabel('y')
        i_row = i_row + 1;
        
        if opt.do_plot_x && ~model.blackbox
            subplot(nrows,1,i_row)
            if strcmp(model.advance_type,'ANN')
                plot(tt,xx(nY+1:end,:),'-','linewidth',1)
            else
                plot(tt,xx,'-','linewidth',1)
            end
            ylabel('x')
        end
        
        pause(1e-16)
    end
end
