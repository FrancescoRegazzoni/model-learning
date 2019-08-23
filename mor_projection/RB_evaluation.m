function ret = RB_evaluation(HFmod,snapshots_source,dataset,opt)
    % Utility to evaluate the convergence of the RB method for the HF model
    % 'HFmod'. It performs POD on the snapshots matrix associated with the
    % source 'snapshots_source', and tests (for different sized of the
    % reduced model, customizable by opt.N_values) the error and the 
    % computational time on the dataset 'dataset'. Optionally, it compares 
    % the results with the models defined in opt.comparison_models_... . It
    % returns the results in a structured object.

    opt.dummy = 0;

    if ~isfield(opt,'N_values')
        opt.N_values = 1:HFmod.nX;
    end
    if ~isfield(opt,'normalize_time')
        opt.normalize_time = 0;
    end    
    if ~isfield(opt,'compute_time_HF')
        opt.compute_time_HF = 0;
    end    
    if ~isfield(opt,'show_epsilon')
        opt.show_epsilon = 0;
    end      
    if ~isfield(opt,'comparison_models_times_nskip')
        opt.comparison_models_times_nskip = 0;
    end
    if ~isfield(opt,'comparison_models_times_ncompute')
        opt.comparison_models_times_ncompute = 1;
    end
    
    if ischar(dataset)
        dataset_def.problem = HFmod.problem;
        dataset_def.type = 'file';
        dataset_def.source = dataset;
        dataset = dataset_get(dataset_def); 
    end
    
    if opt.normalize_time
        if opt.compute_time_HF
            fprintf('Solving HF model to get the mean time...')
            if ischar(opt.dataset_time_HF)                
                dataset_def.problem = HFmod.problem;
                dataset_def.type = 'file';
                dataset_def.source = opt.dataset_time_HF;
                opt.dataset_time_HF = dataset_get(dataset_def);         
            end
            outErrHF = model_compute_error(HFmod,opt.dataset_time_HF);
            if outErrHF.err_dataset_L2_norm > 1e-8
                error('The HF model seems to be wrong (error of %1.2e)...',outErrHF.err_dataset_L2_norm)
            end
            timeRef = outErrHF.time_mean;
        else
            timeRef = opt.mean_time_HF;
        end
        fprintf('NB: time will be normalized by %1.5e s\n',timeRef)
    else
        timeRef = 1;
    end
     
    do_comparison = 0;
    if isfield(opt,'comparison_models')
        if ~isfield(opt.comparison_models{1},'models')
            comparison_models{1}.models = opt.comparison_models;
            comparison_model_names{1} = opt.comparison_model_name;
            opt.comparison_models = comparison_models;
            opt.comparison_model_name = comparison_model_names;
        end
        do_comparison = 1;
        fprintf('Solving comparison models...\n')
        comp_nX=[];
        comp_err=[];
        comp_tim=[];
        for iM = 1:length(opt.comparison_models)
            for iM2 = 1:length(opt.comparison_models{iM}.models)
                time_tmp = 0;
                for i_comp = 1:(opt.comparison_models_times_nskip+opt.comparison_models_times_ncompute)
                    outErr_comparison = model_compute_error(opt.comparison_models{iM}.models{iM2}, dataset);
                    if i_comp > opt.comparison_models_times_nskip
                        time_tmp = time_tmp + outErr_comparison.time_mean;
                    end
                end
                compmod{iM}.comp_nX(iM2) = opt.comparison_models{iM}.models{iM2}.nX;
                compmod{iM}.comp_err(iM2) = outErr_comparison.err_dataset_L2_norm;
                compmod{iM}.comp_tim(iM2) = time_tmp / opt.comparison_models_times_ncompute;
            end
        end
    end  

    %% Build snapshots matrix
    X = build_snapshots_matrix(HFmod.problem,snapshots_source);
    
    %% POD
    
    optPOD.get_full_V = 1;
    optPOD.do_plot = 1;
    optPOD.show_epsilon = 1;
    [Vfull,outputPOD] = POD_projection(X,optPOD);
    
    %% building legend
    fig_err = figure();
    fig_tim = figure();
    fig_errtim = figure();
%     iLeg = 1;
%     legendNames{iLeg} = 'error_{RB}'; iLeg = iLeg + 1;
%     if opt.normalize_time
%         legendNames{iLeg} = 'time_{RB} / time_{HF}'; iLeg = iLeg + 1;
%     else
%         legendNames{iLeg} = 'time_{RB} (s)'; iLeg = iLeg + 1;
%     end
%     if opt.show_epsilon
%         legendNames{iLeg} = '\epsilon'; iLeg = iLeg + 1;
%     end        
%     if do_comparison 
%         legendNames{iLeg} = sprintf('error_{%s}',opt.comparison_model_name); iLeg = iLeg + 1;
%         if opt.normalize_time
%             legendNames{iLeg} = sprintf('time_{%s} / time_{HF}',opt.comparison_model_name); iLeg = iLeg + 1;
%         else
%             legendNames{iLeg} = sprintf('time_{%s} (s)',opt.comparison_model_name); iLeg = iLeg + 1;
%         end
%     end  
    
    iLeg_err = 1;
    iLeg_tim = 1;
    legendNames_err{iLeg_err} = 'error_{RB}'; iLeg_err = iLeg_err + 1;
    if opt.normalize_time
        legendNames_tim{iLeg_tim} = 'time_{RB} / time_{HF}'; iLeg_tim = iLeg_tim + 1;
    else
        legendNames_tim{iLeg_tim} = 'time_{RB} (s)'; iLeg_tim = iLeg_tim + 1;
    end
    if opt.show_epsilon
        legendNames_err{iLeg_err} = '\epsilon'; iLeg_err = iLeg_err + 1;
    end        
    if do_comparison 
        for iM = 1:length(opt.comparison_models)
            legendNames_err{iLeg_err} = sprintf('error_{%s}',opt.comparison_model_name{iM}); iLeg_err = iLeg_err + 1;
            if opt.normalize_time
                legendNames_tim{iLeg_tim} = sprintf('time_{%s} / time_{HF}',opt.comparison_model_name{iM}); iLeg_tim = iLeg_tim + 1;
            else
                legendNames_tim{iLeg_tim} = sprintf('time_{%s} (s)',opt.comparison_model_name{iM}); iLeg_tim = iLeg_tim + 1;
            end
        end
    end  
    
    %% cycle over N
    errs = [];    
    tims = [];      
    epss = [];  
        
    for iN = 1:length(opt.N_values)
        %% building V
        N = opt.N_values(iN);
        V = Vfull(:,1:N);
        epss = [epss POD_getepsilon(outputPOD.s,N)];

        %% Model projection
        RBmod = model_project(HFmod,V,V);

%         %% Analyze the affine dependence of parameters + linearity of equation
%         RBmod_smart = model_makesmart_linear_affineu(RBmod);

        %% Errors
        outErr = model_compute_error(RBmod,dataset);
        errs = [errs outErr.err_dataset_L2_norm];
        tims = [tims outErr.time_mean];
        
        %% plot
%         subplot(2,2,1)
%         plot(opt.N_values(1:iN),errs,'o-')
%         grid on
%         xlabel('N')
%         ylabel('err')
%         
%         subplot(2,2,2)
%         semilogy(opt.N_values(1:iN),errs,'o-')
%         grid on
%         xlabel('N')
%         ylabel('err')
%         
%         subplot(2,2,3)
%         plot(opt.N_values(1:iN),tims,'o-')
%         grid on
%         xlabel('N')
%         ylabel('mean time')

        figure(fig_err)
        semilogy(opt.N_values(1:iN),errs,'o-'); hold on   
        if opt.show_epsilon
            semilogy(opt.N_values(1:iN),epss,'o-')
        end
        
        figure(fig_tim)
        semilogy(opt.N_values(1:iN),tims/timeRef,'o-'); hold on   
        
        figure(fig_errtim)
        semilogy(tims/timeRef,errs,'o-'); hold on   
        
        if do_comparison          
%             ax = gca;
%             ax.ColorOrderIndex = 1;
%             plot(opt.comparison_model.nX,outErr_comparison.err_dataset_L2_norm,'*')
%             plot(opt.comparison_model.nX,outErr_comparison.time_mean/timeRef,'*')
            for iM = 1:length(opt.comparison_models)
                figure(fig_err)
                plot(compmod{iM}.comp_nX,compmod{iM}.comp_err,'*--')
                figure(fig_tim)
                plot(compmod{iM}.comp_nX,compmod{iM}.comp_tim/timeRef,'*--') 
                figure(fig_errtim)
                semilogy(compmod{iM}.comp_tim/timeRef,compmod{iM}.comp_err,'*--');
            end
        end  
    
        figure(fig_err)
        hold off        
        xlabel('N')
        grid on
        legend(legendNames_err)
        
        figure(fig_tim)
        hold off        
        xlabel('N')
        grid on
        legend(legendNames_tim)
        
        figure(fig_errtim)
        hold off        
        xlabel('computational time (relative wrt HF)')
        ylabel('relative error')
        grid on
        %legend(legendNames_errtim)
        
        pause(1e-16)
        
    end
    
    if do_comparison  
        ret.comp_nX = comp_nX;
    end
    ret.compmod = compmod;
    ret.timeRef = timeRef;
    ret.errs = errs;
    ret.tims = tims; 
    if opt.show_epsilon
        ret.epss = epss;
    end
    
end