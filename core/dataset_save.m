function dataset_save(problem, dataset, file_name)

    % Stores a previously generated dataset.

    baseopt = get_base_options();
    full_file_name = [baseopt.BaseDir '/' problem.dir_data '/' file_name];
    
    Tests = dataset;    
    save(full_file_name, 'Tests');
    
end