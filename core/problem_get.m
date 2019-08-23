function problem = problem_get(example,option_file)
    % Gets a problem structure, defined in the file 'option_file' from the
    % folder examples/<example>.

    fprintf('Problem %s - %s loading ...\n',example,option_file);
    
    %% moving into example directory
    currPath = pwd(); 
    baseopt = get_base_options();   
    exampledir = [baseopt.ExampleDir '/' example];
    cd(exampledir);
    problem.example = example;
    problem.goto_dir = @() cd(exampledir);
        
    %% reading options from file    
    problem.name = iniread(option_file,'Problem','name','s','');
    if isempty(problem.name)
         [~,problem.name,~] = fileparts(option_file);
    end
    problem.nU = iniread(option_file,'Problem','nU','d',-1);
    problem.nY = iniread(option_file,'Problem','nY','d',-1);
    problem.y_min = iniread(option_file,'Problem','y_min','d',-1);
    problem.y_max = iniread(option_file,'Problem','y_max','d',1);
    problem.u_min = iniread(option_file,'Problem','u_min','d',-1);
    problem.u_max = iniread(option_file,'Problem','u_max','d',1);
    problem.y0_norm_computation = iniread(option_file,'Problem','y0_norm_computation','d',0);
    
    problem.fixed_x0 = iniread(option_file,'Problem','fixed_x0','d',1);
    problem.samples_variability = iniread(option_file,'Problem','samples_variability','d',0);
    
    % model handlers
    HFmodels = {};
    HFmod_handler_list = iniread(option_file,'Problem','HFmod_handler_list','s','');
    if ~isempty(HFmod_handler_list)
        HFmodels = strsplit(HFmod_handler_list,';');
        for iModel = 1:length(HFmodels)
            CurrentSection = ['Model_' strtrim(HFmodels{iModel})];
            HFmod_handler = iniread(option_file,CurrentSection,'HFmod_handler','s','');
            HFmodels_handlers{iModel}.get_model = eval(HFmod_handler);
        end
        problem.models = containers.Map(HFmodels,HFmodels_handlers);
    end    
    HFmod_handler = iniread(option_file,'Problem','HFmod_handler','s','');
    if isempty(HFmod_handler)
        problem.get_model = @() error('No default model available for this problem');
    elseif any(strcmp(HFmodels,HFmod_handler))
        problem.get_model = problem.models(HFmod_handler).get_model;
    else
        problem.get_model = eval(HFmod_handler);
    end
        
    add_option_if_found('Problem','u_ref','d',-1);
    add_option_if_found('Problem','y_ref','d',-1);
    add_option_if_found('Problem','y0_min','d',-1);
    add_option_if_found('Problem','y0_max','d',-1);
    add_option_if_found('Problem','T','d',-1);
    add_option_if_found('Problem','dir_data','s','');
    add_option_if_found('Problem','dir_nets','s','');
    
    vals = iniread(option_file,'Problem','u_names','s','');
    if ~isequal(vals,'')
        problem.u_names = split(vals,'|');
    end
    vals = iniread(option_file,'Problem','y_names','s','');
    if ~isequal(vals,'')
        problem.y_names = split(vals,'|');
    end
    
    %% launching script if present    
    scriptname = iniread(option_file,'Problem','script','s','');
    if ~isempty(scriptname)
       eval(['problem = ' scriptname ';']);
    end
        
    %% options check
    if problem.nU < 0
        error('Specify a valid value for nU')
    elseif problem.nY < 0
        error('Specify a valid value for nY')
    end
    
    problem = adapt_dimension_struct_field(problem,'y_min'              ,problem.nY);
    problem = adapt_dimension_struct_field(problem,'y_max'              ,problem.nY);
    problem = adapt_dimension_struct_field(problem,'y0_norm_computation',problem.nY);
    problem = adapt_dimension_struct_field(problem,'y_ref'              ,problem.nY);
    problem = adapt_dimension_struct_field(problem,'y0_min'             ,problem.nY);
    problem = adapt_dimension_struct_field(problem,'y0_max'             ,problem.nY);
    
    problem = adapt_dimension_struct_field(problem,'u_min'             ,problem.nU);
    problem = adapt_dimension_struct_field(problem,'u_max'             ,problem.nU);
    problem = adapt_dimension_struct_field(problem,'u_ref'             ,problem.nU);
    
    %% setting authomatic options
    if ~problem.fixed_x0 || problem.samples_variability
        problem.metaproblem = 1; 
        problem.particularized = 0; 
    else
        problem.metaproblem = 0; 
        problem.particularized = 0; 
    end

    %% setting paths
    if ~isfield(problem,'dir_data')
        problem.dir_data = [example '/data_' problem.name ];
    end
    if ~isfield(problem,'dir_nets')
        problem.dir_nets = [example '/networks_' problem.name ];
    end
    create_directory_if_not_found([baseopt.BaseDir '/' example]);
    create_directory_if_not_found([baseopt.BaseDir '/' problem.dir_data]);
    create_directory_if_not_found([baseopt.BaseDir '/' problem.dir_nets]);

    %% moving back to previous directory
    cd(currPath); 
    
    fprintf('Problem %s - %s loaded!\n',example,option_file);
    
    function add_option_if_found(section,option_name,type,null_value)
        val = iniread(option_file,section,option_name,type,null_value);
        if ~isequal(val,null_value)
            eval(['problem.' option_name ' = val;']);
        end
    end

end