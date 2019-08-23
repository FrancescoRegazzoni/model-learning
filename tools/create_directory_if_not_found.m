function create_directory_if_not_found(dirpath)
    % Creates a directory if not existent.

    if ~exist(dirpath,'dir')
        if ~mkdir(dirpath)
            warning('Unable to create directory: %s',dirpath)
        end
    end

end