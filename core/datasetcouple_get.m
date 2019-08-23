function datasetcouple = datasetcouple_get(datasetcouple_def)
    % Gets a training/validation couple of datasets.

    switch datasetcouple_def.type
        case 'file'
            dataset_def.type = 'file';
            dataset_def.problem = datasetcouple_def.problem;
            dataset_def.source = datasetcouple_def.source_train;
            datasetcouple.train = dataset_get(dataset_def);
            dataset_def.source = datasetcouple_def.source_tests;
            datasetcouple.tests = dataset_get(dataset_def);            
        case 'script'
            datasetcouple = GetTrainSet(datasetcouple_def.source);
    end

end