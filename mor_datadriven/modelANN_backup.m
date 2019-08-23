function modelANN_backup(problem, baseDirDest, NetworkDir, idxNet, opt)
    % Backups the model 'NetworkDir' into the folder 'baseDirDest'.

    if nargin == 3
        idxNet = Inf;
    end

    opt.dummy = 0;
    if ~isfield(opt,'name')
        opt.name = NetworkDir;
    end
    if ~isfield(opt,'full')
        opt.full = 1;
    end

    baseopt = get_base_options();

    create_directory_if_not_found([baseDirDest '/' problem.example]);
    create_directory_if_not_found([baseDirDest '/' problem.dir_nets]);    
    create_directory_if_not_found([baseDirDest '/' problem.dir_nets '/' opt.name]);

    pathBase = [baseopt.BaseDir '/' problem.dir_nets];
    pathDest = [baseDirDest     '/' problem.dir_nets];

    iFile = 1;
    filesToCopy{iFile} = 'net.mat'; iFile = iFile+1;
    filesToCopy{iFile} = 'data.mat'; iFile = iFile+1;
    filesToCopy{iFile} = 'normalization.mat'; iFile = iFile+1;
    if opt.full
        filesToCopy{iFile} = 'diary.txt'; iFile = iFile+1;
        fileslist = dir([pathBase '/' NetworkDir '/*.ini']);
        for i = 1:length(fileslist)
            filesToCopy{iFile} = fileslist(i).name; iFile = iFile+1; 
        end
        fileslist = dir([pathBase '/' NetworkDir '/*.m']);
        for i = 1:length(fileslist)
            filesToCopy{iFile} = fileslist(i).name; iFile = iFile+1; 
        end
    else
        filesToCopy{iFile} = 'options.ini'; iFile = iFile+1;
    end

    for i = 1:length(filesToCopy)
        try
            copyfile([pathBase '/' NetworkDir '/' filesToCopy{i}], ...
                     [pathDest '/' opt.name   '/' filesToCopy{i}]);
        catch
            if i <= 2
                error(['unable to copy file ' filesToCopy{i}])
            end
        end    
    end
    if ~isinf(idxNet)
        create_directory_if_not_found([pathDest '/' opt.name '/networks']);
        fileToCopy = sprintf('net_%06d.mat',idxNet);
        copyfile([pathBase '/' NetworkDir '/networks/' fileToCopy], ...
                 [pathDest '/' opt.name   '/networks/' fileToCopy]);

    % 	copyfile([pathBase '/' NetworkDir '/' fileToCopy], ...
    % 			 [pathDest '/' NetworkDir '/networks/' fileToCopy]);
    end

    fprintf('Succesfully copied model %s (idx %d)\n',NetworkDir,idxNet)
    
end