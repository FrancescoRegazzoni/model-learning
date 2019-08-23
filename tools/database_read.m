function db = database_read(filename,problem)
    % Reads a database of trained models.
    
    fprintf('reading database %s...',filename)
    
    db = readtable(filename,'delimiter','|','CommentStyle','%');
    db.problem = strtrim(db.problem);
    db.dir = strtrim(db.dir);
    db.problem = categorical(db.problem);
    if nargin == 2
        db = db(db.problem == problem.name,:);
    end
    
    fprintf(' read!\n')
    
end