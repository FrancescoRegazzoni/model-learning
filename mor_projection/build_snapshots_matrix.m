function X = build_snapshots_matrix(problem, source)
    % Builds the snapshots matrix for the problem 'problem', starting from
    % the states stored into the tests associated with 'source'.

    fprintf('building snapshots matrix...')
    X = [];
    baseopt = get_base_options();
    
    sources = strsplit(source,'|');
    for i = 1:length(sources)
        fields = strsplit(sources{i},';');
        filename = [baseopt.BaseDir '/' problem.dir_data '/' fields{1}];
        if exist(filename, 'file') == 2
            SamplesGen = load(filename);
            if strcmp(fields{2},':')
                idxs = 1:length(SamplesGen.Tests);
            else
                idxs = eval(fields{2});
            end
            for isampl = 1:length(idxs)
                X = [X SamplesGen.Tests{idxs(isampl)}.xx];
            end
            clear SamplesGen
        else
            error('file %s does not exist!',filename)
        end
    end
    
    fprintf(' done!\n')
end