function Tests = dataset_generate(model,dataset,opt)

    % Generates a collection of tests, by solving 'model' on 'dataset'.
    % Optionally, it stores the results in the file opt.outFile.

    %% options
    opt.dummy = 0;
    
    if ~isfield(opt,'do_save')
        opt.do_save = 0;
    end
    if ~isfield(opt,'outFile')
        opt.outFile = '';
    end
    if ~isfield(opt,'freq_save')
        opt.freq_save = 10;
    end
    if ~isfield(opt,'do_plot')
        opt.do_plot = 0;
    end
    
    %% initialization
    if opt.do_save 
        baseopt = get_base_options();
        filename = [baseopt.BaseDir '/' model.problem.dir_data '/' opt.outFile];
    end
    
    %% dataset generation
    for iTest = 1:length(dataset)
        fprintf('solving test %d of %d...\n',iTest,length(dataset))
        Tests{iTest} = model_solve(dataset{iTest},model,struct('verbose',1));  
        if opt.do_plot
            if iTest > 0, close; end
            dataset_plot({Tests{iTest}},model.problem);
            pause(1e-16)
        end
        
        %% saving file
        if opt.do_save && ( mod(iTest,opt.freq_save)==0 || iTest == length(dataset) )
            save(filename,'Tests');
        end
    end
end