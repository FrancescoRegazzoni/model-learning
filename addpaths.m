function addpaths(mainpath)
    % Adds library folders to Matlab paths.

    if nargin == 0
        [mainpath,~,~] = fileparts(mfilename('fullpath'));
    end

    addpath([mainpath '/core'])
    addpath([mainpath '/tools'])
    addpath([mainpath '/ANN'])
    addpath([mainpath '/optimization'])
    addpath([mainpath '/mor_projection'])
    addpath([mainpath '/mor_datadriven'])
    addpath([mainpath '/data_assimilation'])

    % addpath(genpath(...));
    % rmpath(genpath(...));

    fprintf(' --- model learning: path added.\n')

end