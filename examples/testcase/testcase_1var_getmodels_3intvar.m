function model = testcase_1var_getmodels_3intvar(problem, name)

    fprintf('building toy problem model...')
    model.problem = problem;

    model.nA = 1;
    model.alpha_min = 0;
    model.alpha_max = 1;

    model.nX = 3;
    model.x0_min = [1  ;-.05;-.05];
    model.x0_max = [1.1; .05; .05];
%     model.x0_min = [1  ;1;1];
%     model.x0_max = [1.1;1;1];
    model.dt = 1e-4;
    model.advance_type = 'nonlinear_explicit';
    model.u_implicit = 0;
    model.output_type = 'nonlinear';
	model.g = @(x) x(1)+x(2);
	model.dgdx = @(x) [1,1,0]; 
	
	c = 0;
	period_osc = .1;
	additionaldynamics = @(x) [x(2); -c*x(2)-(2*pi/period_osc)^2*x(1)];
	dadditionaldynamics = @(x) [0, 1; -(2*pi/period_osc)^2, -c];
    
    switch name
        case 'linear'
            f_alpha = @(x,u,a) 2*a;
            dfdx = @(x,u,a) 0;
            dfda = @(x,u,a) 2;
        case 'exp'
            f_alpha = @(x,u,a) a*x;
            dfdx = @(x,u,a) a;
            dfda = @(x,u,a) x;
        case 'expstab'
            f_alpha = @(x,u,a) 2*a*(3 - x);
            dfdx = @(x,u,a) -2*a;
            dfda = @(x,u,a) 2*(3 - x);
        case 'bistab'
            f_alpha = @(x,u,a) 2*a*x*(3 - x);
            dfdx = @(x,u,a) 6*a - a*x;
            dfda = @(x,u,a) 2*x*(3 - x);
    end
    
    model.f_alpha = @(x,u,a) [f_alpha(x(1),u,a); additionaldynamics(x(2:3))];
    model.dfdx = @(x,u,a) [dfdx(x(1),u,a), 0, 0; [0;0], dadditionaldynamics(x(2:3))];
    model.dfda = @(x,u,a) [dfda(x(1),u,a);0;0];
    
    fprintf(' done!\n')
    
%     model.particularize_model = @particularize_model;
%     function mod_part = particularize_model(mod,alpha)
%         mod_part = mod;
%         mod_part.alpha = alpha;
%         
%         switch name
%             case 'linear'
%                 mod_part.f = @(x,u) [2*mod_part.alpha;                  additionaldynamics(x(2:3))];
%             case 'exp'
%                 mod_part.f = @(x,u) [mod_part.alpha*x(1);               additionaldynamics(x(2:3))];
%             case 'expstab'
%                 mod_part.f = @(x,u) [2*mod_part.alpha*(3 - x(1));       additionaldynamics(x(2:3))];
%             case 'bistab'
%                 mod_part.f = @(x,u) [2*mod_part.alpha*x(1)*(3 - x(1));  additionaldynamics(x(2:3))];
%         end
%     
% %         mod_part.f = @(x,u) [0;additionaldynamics(x(2:3))];
% %         mod_part.f = @(x,u) [mod_part.alpha*x(1);x(3);-c*x(3)-(2*pi/period_osc)^2*x(2)];
%         
%     end
    
end