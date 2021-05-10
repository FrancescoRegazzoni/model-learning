function st = steady_state_load(problem,filename)
    % Loads the steady-states stores into the file 'filename'.
    
    baseopt = get_base_options();
    filename = [baseopt.BaseDir '/' problem.dir_data '/' filename];
    fprintf('loading dataset... ');
    st = load(filename);
    fprintf('done!\n');
    
end