function dataset_export_table(problem, dataset, foldername)

    create_directory_if_not_found(foldername)
    
    for iS = 1:length(dataset)
        filename_u = sprintf('%s/sample_%d_input.dat', foldername, iS);
        filename_y = sprintf('%s/sample_%d_output.dat', foldername, iS);
        if isfield(dataset{iS}, 'uu')
            lables = [{'t'}, problem_getvariablename_list(problem, 'u')];
            table = array2table([dataset{iS}.tt' dataset{iS}.uu'],lables);
            writetable(table,filename_u,'Delimiter','tab');
        end
        
        if ~isfield(dataset{iS}, 'tt_y')
            dataset{iS}.tt_y = dataset{iS}.tt;
        end
        lables = [{'t'}, problem_getvariablename_list(problem, 'y')];
        table = array2table([dataset{iS}.tt_y' dataset{iS}.yy'],'VariableNames',lables);
        writetable(table,filename_y,'Delimiter','tab')
        
        fprintf('exported sample %d of %d\n', iS, length(dataset))
    end
end