function da_test_on_training_set(ANNmod, obs_err, mod_err)

    dataset_train = getfield(datasetcouple_get(ANNmod.datasets_def), 'train');
    figure('units','normalized','outerposition',[0 0 1 1])
    for i = 1:length(dataset_train)
        out_da = model_da(ANNmod,dataset_train{i},struct('obs_err',obs_err,'mod_err',mod_err));
        da_opt.dummy = 0;
        if ANNmod.problem.samples_variability
            da_opt.alpha_true = ANNmod.alpha_learned(:,i);
        end
        da_plot_output(out_da,dataset_train{i},ANNmod,da_opt);
        pause()
    end

end