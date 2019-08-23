function out = iniread(file,section,field,type,defaultval)
    % Reads the content of an ini file (type = 'i','d','s')
    
    if nargin == 3
        type = 'd';
    end
    key = {section,'',field,type};
    if nargin == 5
        key{5} = defaultval;
    end
    sTemp = inifile(file,'read',key); 
    out = sTemp{1};
end