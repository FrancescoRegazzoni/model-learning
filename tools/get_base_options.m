function opt = get_base_options()
    % Returns an object containing the base options of the library.
    
    currPath = pwd();
    [filepath,~,~] = fileparts(mfilename('fullpath'));
    [codepath,~,~] = fileparts(filepath);

    cd(codepath);

    if isfile('options.ini')
        opt.BaseDir = iniread('options.ini','paths','datapath','s');
    else
        opt.BaseDir = '%/data';
    end
    if opt.BaseDir(1) == '%'
        opt.BaseDir = [codepath opt.BaseDir(2:end)];
    end
    opt.CodeDir = codepath;
    opt.ExampleDir = [codepath '/examples'];

    cd(currPath);

end