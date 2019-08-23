function name = problem_getvariablename(problem,type,index)

    switch type
        case {'u','y'}
            if isfield(problem,[type '_names'])
                names = getfield(problem,[type '_names']);
                name = names{index};
            else
                if (isequal(type,'u') && problem.nU == 1) || (isequal(type,'y') && problem.nY == 1)
                    name = type;
                else
                    name = sprintf([type '_%d'],index);
                end
            end
        otherwise
            error('variable type not recognized')
    end

end