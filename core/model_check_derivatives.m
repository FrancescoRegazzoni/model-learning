function model_check_derivatives(mod)
    % Check the derivatives associated with the definition of a model.

    ep_rel = 1e-6;
    ep_min = 1e-8;
    n_samples = 1e3;
    
    if isfield(mod,'f_alpha')
        check_dfdx = isfield(mod,'dfdx');
        check_dfda = isfield(mod,'dfda');
        if check_dfdx, errtot_dfdx = 0; end
        if check_dfda, errtot_dfda = 0; end
        for i = 1:n_samples
            mod_part = metamodel_particularize(mod);
            x = mod_part.x0;
            a = mod_part.alpha;
            if mod.problem.nU == 0
                u = [];
            else
                u = mod.problem.u_min + (mod.problem.u_max-mod.problem.u_min)*rand(1);
            end
            f_0 = mod.f_alpha(x,u,a);
            if check_dfdx
                dfdx_ana = mod.dfdx(x,u,a);
                dfdx_num = zeros(mod.nX,mod.nX);
                for iX = 1:mod.nX
                    ep = max(ep_min,abs(x(iX))*ep_rel);
                    x_ep = x;
                    x_ep(iX) = x_ep(iX) + ep;
                    f_ep = mod.f_alpha(x_ep,u,a);
                    dfdx_num(:,iX) = (f_ep-f_0)/ep;
                end
                errtot_dfdx = errtot_dfdx + norm(dfdx_num-dfdx_ana)/max([norm(dfdx_num),norm(dfdx_ana),eps]);
            end
            if check_dfda
                dfda_ana = mod.dfda(x,u,a);
                dfda_num = zeros(mod.nX,mod.nA);
                for iA = 1:mod.nA
                    ep = max(ep_min,abs(a(iA))*ep_rel);
                    a_ep = a;
                    a_ep(iA) = a_ep(iA) + ep;
                    f_ep = mod.f_alpha(x,u,a_ep);
                    dfda_num(:,iA) = (f_ep-f_0)/ep;
                end
                errtot_dfda = errtot_dfda + norm(dfda_num-dfda_ana)/max([norm(dfda_num),norm(dfda_ana),eps]);
            end
        end
        fprintf('num samples: %d, ep_rel = %1.0e, ep_min = %1.0e\n',n_samples,ep_rel,ep_min)
        if check_dfdx
            fprintf('dfdx rel err: %1.2e\n',errtot_dfdx/n_samples)
        end
        if check_dfda
            fprintf('dfda rel err: %1.2e\n',errtot_dfda/n_samples)
        end
    end
end