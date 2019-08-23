function model_smart = model_makesmart_linear_affineu(model)
    % A linear model, can be written in the form:
    %     K x^(k+1) = B x^(k) + b
    % where K,B,b depend on dt (which we suppose to be given)
    % Moreovoer, in case of affine input dependence, we have:
    %     K = K_0 + u(1)*K_1 + u(2)*K_2 + ...
    %     B = B_0 + u(1)*B_1 + u(2)*B_2 + ...
    %     b = b_0 + u(1)*b_1 + u(2)*b_2 + ...
    % By calling the "solveonestep_linear" function of the original model
    % several times, each time with a different element of a basis of both
    % x and u, it is possible to reconstruct the K_j, B_j and b_j.
    
    fprintf('Affine dependence analysis:\n')
    
    if ~strcmp(model.advance_type,'solveonestep_linear')
        error('This utility is only for models of type ''solveonestep_linear''!')
    end
    
    model_smart = model;
    
    nU = model.problem.nU;
    nX = model.nX;
    dt = model.dt;
    
    [model_smart.K0,model_smart.b0,model_smart.B0] = get_KbB(0);
    model_smart.K = zeros(nX,nX,nU);
    model_smart.b = zeros(nX,nU);
    model_smart.B = zeros(nX,nX,nU);
    for iU = 1:nU
        fprintf('   analyzing parameter %d/%d...\n',iU,nU)
        [K,b,B] = get_KbB(iU);
        model_smart.K{iU} = K - model_smart.K0;
        model_smart.b{iU} = b - model_smart.b0;
        model_smart.B{iU} = B - model_smart.B0;
%         model_smart.K(:,:,iU) = K - model_smart.K0;
%         model_smart.b(:,iU) =   b - model_smart.b0;
%         model_smart.B(:,:,iU) = B - model_smart.B0;
    end
    
%     model_smart.solveonestep_linear = @solveonestep_linear_smart;
    model_smart.advance_type = 'linear_advance_uaffine';
    
    fprintf('done!\n')
    
    function [K,b,B] = get_KbB(idxU)
        
        u = zeros(nU,1);
        if idxU>0
            u(idxU) = 1;
        end
        
        x = zeros(nX,1);        
        [K,b] = model.solveonestep_linear(x,u,dt);
        
        B = zeros(nX,nX);
        for iX = 1:nX
            x = zeros(nX,1);
            x(iX) = 1;
            [~,rhs] = model.solveonestep_linear(x,u,dt);
            B(:,iX) = rhs - b;
        end
    end

%     function [K,rhs] = solveonestep_linear_smart(x,u,dt)
%         if nargin == 3
%             if abs(dt - model.dt)/dt > 1e-3
%                 error('using time step different than that used to make the model smart')
%             end
%         end
%         K = model_smart.K0;
%         rhs = model_smart.b0 + model_smart.B0*x;
%         for iU = 1:nU
%             K = K + u(iU)*model_smart.K(:,:,iU);
%             rhs = rhs + u(iU)*(model_smart.b(:,iU) + model_smart.B(:,:,iU)*x);
%         end
%     end

end