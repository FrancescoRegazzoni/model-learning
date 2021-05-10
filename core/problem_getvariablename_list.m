function names = problem_getvariablename_list(problem,type)

    if nargin == 1
        type = 'a';
    end
    switch type
        case {'u','y'}
            if isequal(type,'u'), nvar = problem.nU; end
            if isequal(type,'y'), nvar = problem.nY; end
            if nvar == 0
                names = [];
            else
                for i = 1:nvar
                    names{i} = problem_getvariablename(problem,type,i);
                end
            end
        case 'a'
            names = [ problem_getvariablename_list(problem,'u'), ...
                      problem_getvariablename_list(problem,'y') ];
        otherwise
            error('variable type not recognized')
    end

end