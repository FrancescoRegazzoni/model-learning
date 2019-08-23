function mod_part = metamodel_particularize(mod,x0,alpha)
    % Particularizes a meta-model, wrt initial condition (x0) and/or model 
    % parameters (alpha). If not needed, or to use random initialization, 
    % don't pass the input or use an empty vector.

    if nargin < 2
        x0 = [];
    end
    if nargin < 3
        alpha = [];
    end
    
    mod_part = mod;
    mod_part.problem.particularized = 1;
    
    if ~mod.problem.fixed_x0
        if isempty(x0)
            if isfield(mod,'x0_getrandom')
                mod_part.x0 = mod.x0_getrandom();
            else
                mod_part.x0 = mod.x0_min + (mod.x0_max - mod.x0_min).*rand(mod.nX,1);
            end
        else
            mod_part.x0 = x0;
        end
    end
    
    if mod.problem.samples_variability
        if isempty(alpha)
            if isfield(mod,'alpha_getrandom')
                mod_part.alpha = mod.alpha_getrandom();
            else
                mod_part.alpha = mod.alpha_min + (mod.alpha_max - mod.alpha_min).*rand(mod.nA,1);
            end
        else
            mod_part.alpha = alpha;
        end
        
        if isfield(mod,'particularize_model')
            mod_part = mod_part.particularize_model(mod_part,mod_part.alpha);
        else
            switch mod.advance_type
                case 'nonlinear_explicit'
                    mod_part.f = @(x,u) mod_part.f_alpha(x,u,mod_part.alpha);
                case 'linear_advance_timeexplicit_uaffine'                    
                    if isfield(mod_part,'Bf0_alpha'), mod_part.Bf0 = mod_part.Bf0_alpha(mod_part.alpha); end
                    if isfield(mod_part,'Kf0_alpha'), mod_part.Kf0 = mod_part.Kf0_alpha(mod_part.alpha); end
                    if isfield(mod_part,'Bv0_alpha'), mod_part.Bv0 = mod_part.Bv0_alpha(mod_part.alpha); end
                    if isfield(mod_part,'Kv0_alpha'), mod_part.Kv0 = mod_part.Kv0_alpha(mod_part.alpha); end
                    if isfield(mod_part,'b0_alpha' ), mod_part.b0  = mod_part.b0_alpha( mod_part.alpha); end
                    if isfield(mod_part,'Bf_alpha'), for i=1:mod_part.problem.nU, mod_part.Bf{i} = mod_part.Bf_alpha{i}(mod_part.alpha); end, end
                    if isfield(mod_part,'Kf_alpha'), for i=1:mod_part.problem.nU, mod_part.Kf{i} = mod_part.Kf_alpha{i}(mod_part.alpha); end, end
                    if isfield(mod_part,'Bv_alpha'), for i=1:mod_part.problem.nU, mod_part.Bv{i} = mod_part.Bv_alpha{i}(mod_part.alpha); end, end
                    if isfield(mod_part,'Kv_alpha'), for i=1:mod_part.problem.nU, mod_part.Kv{i} = mod_part.Kv_alpha{i}(mod_part.alpha); end, end
                    if isfield(mod_part,'b_alpha' ), for i=1:mod_part.problem.nU, mod_part.b{i}  = mod_part.b_alpha{ i}(mod_part.alpha); end, end 
                case 'linear_CN_uaffine'
                    mod_part.A = mod_part.A_alpha(mod_part.alpha);
                    if isfield(mod_part,'b0_alpha' ), mod_part.b0  = mod_part.b0_alpha( mod_part.alpha); end
                    if isfield(mod_part,'b_alpha' )
                        for i=1:mod_part.problem.nU
                            mod_part.b{i}  = mod_part.b_alpha{ i}(mod_part.alpha); 
                        end
                    end
                otherwise
                    error('not (yet) implemented')
            end
        end
        
    end    
end