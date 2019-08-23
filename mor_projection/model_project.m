function model_projected = model_project(model,V,W)
    % Computes the projection-based reduced model by means of the matrices
    % V and W.

    fprintf('Projecting model...')
    
    if isfield(model,'get_projected_model')
        model_projected = model.get_projected_model(V,W);
    else
            
        model_projected = model;
        
        %% Projecting the initial condition
        model_projected.x0 = V'*model.x0;
        model_projected.nX = length(model_projected.x0);
        
        %% Projecting the output
        switch model.output_type
            case 'linear' 
                model_projected.G = model.G*V;
            case 'nonlinear'
                model_projected.g = @(x) model.g(V*x);
            otherwise
                error('unknown output type %s',model.output_type)
        end
        
        %% Projecting the equation
        switch model.advance_type
            case 'solveonestep_linear'
                model_projected.solveonestep_linear = @solveonestep_linear_projected;
            case 'nonlinear_explicit'
                model_projected.f = @(x,u) W'*model.f(V*x,u);
            case 'ANN'
                error('not (yet) implemented!')
            case 'linear_advance'
                model_projected.solveonestep_linear = @linear_advance_projected;
            case 'linear_advance_timeexplicit'
                model_projected.solveonestep_linear = @linear_advance_timeexplicit_projected;
            case 'linear_advance_uaffine'
                model_projected.K0 = W'*model.K0*V;
                model_projected.B0 = W'*model.B0*V;
                model_projected.b0 = W'*model.b0;
                for iU = 1:model.problem.nU
                    model_projected.K{iU} = W'*model.K{iU}*V;
                    model_projected.B{iU} = W'*model.B{iU}*V;
                    model_projected.b{iU} = W'*model.b{iU};
    %                 model_projected.K(:,:,iU) = W'*model.K(:,:,iU)*V;
    %                 model_projected.B(:,:,iU) = W'*model.B(:,:,iU)*V;
    %                 model_projected.b(:,iU)   = W'*model.b(:,iU);
                end
            case 'linear_advance_timeexplicit_uaffine'
                model_projected.Kf0 = W'*model.Kf0*V;
                model_projected.Kv0 = W'*model.Kv0*V;
                model_projected.Bf0 = W'*model.Bf0*V;
                model_projected.Bv0 = W'*model.Bv0*V;
                model_projected.b0  = W'*model.b0;
                for iU = 1:model.problem.nU
                    model_projected.Kf{iU} = W'*model.Kf{iU}*V;
                    model_projected.Kv{iU} = W'*model.Kv{iU}*V;
                    model_projected.Bf{iU} = W'*model.Bf{iU}*V;
                    model_projected.Bv{iU} = W'*model.Bv{iU}*V;
                    model_projected.b{iU}  = W'*model.b{iU};
    %                 model_projected.Kf(:,:,iU) = W'*model.Kf(:,:,iU)*V;
    %                 model_projected.Kv(:,:,iU) = W'*model.Kv(:,:,iU)*V;
    %                 model_projected.Bf(:,:,iU) = W'*model.Bf(:,:,iU)*V;
    %                 model_projected.Bv(:,:,iU) = W'*model.Bv(:,:,iU)*V;
    %                 model_projected.b(:,iU)    = W'*model.b(:,iU);
                end
            otherwise
                error('Impossible to project a model of this kind')
        end
        
    end
    
    fprintf(' done!\n')
    
    function [K,rhs] = solveonestep_linear_projected(x,u,dt)
        [K_full,rhs_full] = model.solveonestep_linear(V*x,u,dt);
        K = W'*K_full*V;
        rhs = W'*rhs_full;
    end

    function [K,B,b] = linear_advance_projected(u,dt)
        [K_full,B_full,b_full] = model.linear_advance(u,dt);
        K = W'*K_full*V;
        B = W'*B_full*V;
        b = W'*b_full;
    end

    function [Kf,Kv,Bf,Bv,b] = linear_advance_timeexplicit_projected(u)        
       [Kf_full,Kv_full,Bf_full,Bv_full,b_full] = model.linear_advance_timeexplicit(u);
        Kf = W'*Kf_full*V;
        Kv = W'*Kv_full*V;
        Bf = W'*Bf_full*V;
        Bv = W'*Bv_full*V;
        b = W'*b_full;
    end

end