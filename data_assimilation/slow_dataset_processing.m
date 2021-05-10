function Tests = slow_dataset_processing(pb_slow,ds_slow,mod_fast,opt)

%% Default options
opt.dummy = 0;
if ~isfield(opt,'noise_abs')
    opt.noise_abs = 0;
end
if ~isfield(opt,'noise_rel_wrt_mean')
    opt.noise_rel_wrt_mean = 0;
end
if ~isfield(opt,'noise_rel_wrt_range')
    opt.noise_rel_wrt_range = 0;
end
if ~isfield(opt,'outFile')
    opt.outFile = '';
end
if ~isfield(opt,'da_obs_noise')
    opt.da_obs_noise = -1;
end
if ~isfield(opt,'do_plot')
    opt.do_plot = 0;
end
if ~isfield(opt,'pause_each_test')
    opt.pause_each_test = 0;
end

%% Initialization
if ischar(ds_slow) 
    dataset_def.problem = pb_slow;
    dataset_def.type = 'file';
    dataset_def.source = ds_slow;
    ds_slow = dataset_get(dataset_def);
end

if isfield(ds_slow{1},'yy')
    synthetic_data = 1;
else
    synthetic_data = 0;
end

Tests = ds_slow;

y_mean = (mod_fast.problem.y_max+mod_fast.problem.y_min)/2;
y_range = (mod_fast.problem.y_max-mod_fast.problem.y_min)/2;
noise_magnitude = opt.noise_abs ...
              + opt.noise_rel_wrt_range.*y_range ... 
              + opt.noise_rel_wrt_mean.*y_mean ;
noise_magnitude = adapt_dimension(noise_magnitude,mod_fast.problem.nY);
noise_vals = [noise_magnitude noise_magnitude./y_mean noise_magnitude./y_range];
fprintf('Noise magnitude (abs --- rel wrt ymean --- rel wrt range):\n')
fprintf('%e --- %e --- %e\n',noise_vals)
fprintf('Error magnitude (abs --- rel wrt ymean --- rel wrt range):\n')
fprintf('%e --- %e --- %e\n',noise_vals/sqrt(mod_fast.dt))

if isequal(opt.da_obs_noise,-1)
    da_opt.obs_err = noise_magnitude;
end
da_opt.obs_err = opt.da_obs_noise;
% da_opt.do_plot = opt.do_plot;
if opt.do_plot
    figure()
end

%% Loop over samples and over slow time
for iS = 1:length(ds_slow)
    fprintf('------- sample %d/%d\n',iS,length(ds_slow))
    
    dt_slow = ds_slow{iS}.tt(2:end)-ds_slow{iS}.tt(1:end-1);
    dt_slow = [dt_slow dt_slow(end)];
    
    tic
    for jTau=1:length(ds_slow{iS}.tt)
        fprintf(' time %d/%d ... ',jTau,length(ds_slow{iS}.tt))
        
        if (synthetic_data)            
            %% particularize model
            fast_params_ij = ds_slow{iS}.yy(:,jTau);
            Tests{iS}.yy_ex(:,jTau) = fast_params_ij;
            mod_fast_ij = metamodel_particularize(mod_fast,[],fast_params_ij);
            % NB: in this case we are assuming that the fast params
            % coincide with the slow state. In case the slow model have
            % some further option which can depend on patient&time, only
            % time or constant, they may be stored in the dataset and
            % loaded here
            
            %% solve fast model
            if isfield(ds_slow{iS},'fast_tests') 
                % patient & time specific fast scale input
                fast_test_ij = ds_slow{iS}.fast_tests{jTest};
            elseif isfield(ds_slow{iS},'fast_test_base') 
                % patient specific fast scale input
                fast_test_ij = ds_slow{iS}.fast_test_base;
            else
                % generic fast scale input
                fast_test_ij = opt.fast_test_base;
            end
            fast_test_ij = model_solve(fast_test_ij,mod_fast_ij);
            
            %% add noise
            fast_test_ij.yy = fast_test_ij.yy + noise_magnitude/sqrt(mod_fast.dt).*randn(size(fast_test_ij.yy));
        else
            fast_test_ij = ds_slow{iS}.fast_tests{jTest};
        end
        
        %% data assimilation 
        out_da = model_da(mod_fast,fast_test_ij,da_opt);
        fast_params_ij_hat = out_da.xx(mod_fast.nX+1:end,end);    
        Tests{iS}.yy(:,jTau) = fast_params_ij_hat;
        
        if synthetic_data
            error_abs = sqrt(mean((fast_params_ij-fast_params_ij_hat).^2));
            error_rel = error_abs/sqrt(mean(fast_params_ij.^2));         
            noise_abs = error_abs*sqrt(dt_slow(jTau));        
            noise_rel = error_rel*sqrt(dt_slow(jTau));
            fprintf('error: (abs) %1.4f, (rel) %1.4f --- noise (abs) %1.4f, (rel) %1.4f\n',error_abs,error_rel,noise_abs,noise_rel)
        end
        
        if opt.do_plot
            opt_plot.dummy = 0;
            if synthetic_data
                opt_plot.alpha_true = fast_params_ij;
            end
            da_plot_output(out_da,fast_test_ij,mod_fast,opt_plot)
            if opt.pause_each_test
                pause()
            else
                pause(1e-16)
            end
        end

    end
    toc
    
    error_abs = sqrt(mean((Tests{iS}.yy(:) - Tests{iS}.yy_ex(:)).^2));
    error_rel = error_abs/sqrt(mean(Tests{iS}.yy_ex(:).^2));    
    fprintf('mean values over sample %d:\n',iS)
    fprintf('error: (abs) %1.4f, (rel) %1.4f \n',error_abs,error_rel)
end

if ~isempty(opt.outFile)
    baseopt = get_base_options();
    filename = [baseopt.BaseDir '/' pb_slow.dir_data '/' opt.outFile '.mat'];
    save(filename,'Tests')
end