function dataset_generate_trace(pb_full,ds_full,pb_trace,ds_trace,trace_function_x,trace_function_y,trace_function_a)

    fprintf('Taking the trace of %s...\n',ds_full)
    
    dataset_def.problem = pb_full;
    dataset_def.type = 'file';
    dataset_def.source = ds_full;
    Tests = dataset_get(dataset_def);
    
    for iS = 1:length(Tests)
        Tests{iS}.x0 = trace_function_x(Tests{iS}.x0);
        Tests{iS}.alpha = trace_function_a(Tests{iS}.alpha);
        Tests{iS}.yy = trace_function_y(Tests{iS}.yy);
        if isfield(Tests{iS},'yy_ex')
            Tests{iS}.yy_ex = trace_function_y(Tests{iS}.yy_ex);
        end
    end
    
    baseopt = get_base_options();
    filename = [baseopt.BaseDir '/' pb_trace.dir_data '/' ds_trace '.mat'];
    save(filename,'Tests')
    
    fprintf('Trace saved into %s\n',ds_trace)
    
end