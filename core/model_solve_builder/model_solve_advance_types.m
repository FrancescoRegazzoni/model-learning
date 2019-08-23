%%CASE 'solveonestep' 
	% the most generic one                                                                 
	x = model.solveonestep(x,uu_eff(:,iT),dtt(iT));                                          
%%CASE 'solveonestep_linear'                                                                             
	% like "solveonestep" but this entails the solution of a linear system                 
	[K,rhs] = model.solveonestep_linear(x,uu_eff(:,iT),dtt(iT));                             
	x = K \ rhs;                                                                            
%%CASE 'nonlinear_explicit'                                                                              
	% explicit euler, with generic right-hand side (NB: model.u_implicit should be false!) 
	x = x + dtt(iT) * model.f(x,uu_eff(:,iT));                                               
%%CASE 'ANN'                                                                                             
	% explicit euler, with ANN right-hand side (NB: model.u_implicit should be false!)     
	if nU > 0                                                                                
		x = x + dtt(iT) * model.ANN([uu_eff(:,iT);x]);                                      
	else                                                                                     
		x = x + dtt(iT) * model.ANN(x);                                                     
	end                                                                                      
%%CASE 'linear_CN_uaffine'                                                                               
   % generic linear model, with Crank-Nicolson, with affine dependence on u                
   % dx/dt = A x + (b0 + sum_j b_j u_j)                                                    
   u = uu_eff(:,iT);                 
   b  = model.b0;                                                                           
   for iU = 1:nU                                                                            
	   b  = b  + u(iU)*model.b{iU};                                                         
   end                                                                                      
   x = (eye(nX)/dtt(iT) - .5*model.A) \ ((eye(nX)/dtt(iT) + .5*model.A)*x + b);            
%%CASE 'linear_advance'                                                                                  
	% generic linear model                                                                 
	[K,B,b] = model.linear_advance(uu_eff(:,iT),dtt(iT));                                    
	x = K \ (B*x + b);                                                                      
%%CASE 'linear_advance_timeexplicit'                                                                     
	% generic linear model, with explicit dependence on dt                                 
	[Kf,Kv,Bf,Bv,b] = model.linear_advance_timeexplicit(uu_eff(:,iT));                       
	x = (Kf + Kv/dtt(iT)) \ ((Bf + Bv/dtt(iT))*x + b);                                      
%%CASE 'linear_advance_uaffine'                                                                          
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
%%CASE 'linear_advance_timeexplicit_uaffine'                                                             
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