function backup_data(problem, baseDirDest, datasets)
    % Backups the data 'datasets' into the folder 'baseDirDest'.

    baseopt = get_base_options();

    create_directory_if_not_found([baseDirDest '/' problem.example]);
    create_directory_if_not_found([baseDirDest '/' problem.dir_data]);

    pathBase = [baseopt.BaseDir '/' problem.dir_data];
    pathDest = [baseDirDest     '/' problem.dir_data];

    for i = 1:length(datasets)
        copyfile([pathBase '/' datasets{i}], ...
                 [pathDest '/' datasets{i}]);
             
        fprintf('Succesfully copied dataset %s\n',datasets{i})
    end

end