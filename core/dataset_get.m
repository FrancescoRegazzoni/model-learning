function dataset = dataset_get(dataset_def)
    % Gets a dataset.

    fprintf('loading dataset...')    
%     dataset.problem = dataset_def.problem;
    
    switch dataset_def.type
        case 'file'
            % dataset_def.source = file1;A;B;C|file2;A|... 
            % A = range (default: ":", which stands for "all")
            % B = 0/1 closed loop (optional, default: 0)
            % C = 0/1 closed loop (optional, default: 1)
            % D = dt (optional)
            
            if isempty(dataset_def.source)
                dataset = [];
            else
                sources = strsplit(dataset_def.source,'|');
                icurr = 1;
                baseopt = get_base_options();
                for i = 1:length(sources)
                    fields = strsplit(sources{i},';');
                    filename = [baseopt.BaseDir '/' dataset_def.problem.dir_data '/' fields{1}];
                    if length(fields) == 1
                        fields{2} = ':';
                    end
                    if exist(filename, 'file') == 2
                        SamplesGen = load(filename);
                        if strcmp(fields{2},':')
                            idxs = 1:length(SamplesGen.Tests);
                        else
                            idxs = eval(fields{2});
                        end
                        for isampl = 1:length(idxs)
                            dataset{icurr} = SamplesGen.Tests{idxs(isampl)};

                            if length(fields) >= 3
                                dataset{icurr}.closed_loop = str2double(fields{3});
                            else
                                dataset{icurr}.closed_loop = 0;
                            end
                            if length(fields) >= 4
                                dataset{icurr}.eval_diff = str2double(fields{4});
                            else
                                dataset{icurr}.eval_diff = 1;
                            end
                            if length(fields) >= 5
                                dataset{icurr}.dt_eval = fields{5};
                            end

                            icurr = icurr+1;
                        end
                        clear SamplesGen
                    else
                        error('file %s does not exist!',filename)
                    end
                end
            end

        case 'script'

            trainset = GetTrainSet(dataset_def.source);
            if strcmp(dataset_def.which_set,'train')
                dataset = trainset.train;
            else
                dataset = trainset.tests;
            end

    end

    fprintf(' done!\n')
    
end