function modelANN_backupmodels(ANNmodels, baseDirDest, opt)
    % Backups the models 'ANNmodels' into the folder 'baseDirDest'.
    
    opt.dummy = 0;
    
    for imod = 1:length(ANNmodels)
        if isfield(opt,'names')
            opt.name = opt.names{imod};
        end
        modelANN_backup(ANNmodels{imod}.problem, baseDirDest, ANNmodels{imod}.NetworkDir, ANNmodels{imod}.idxNet, opt);
    end
end