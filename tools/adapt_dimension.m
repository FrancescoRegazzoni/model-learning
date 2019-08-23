function output = adapt_dimension(input,dim)
    % Transforms 'intput' into a dim-by-1 matrix.
    if dim == 0
        output = [];
    elseif isequal(size(input),[dim,1])
        output = input;
    elseif isequal(size(input),[1,1])
        output = input * ones(dim,1);
    elseif isequal(size(input),[1,dim])
        output = input';
    else
        error('incompatible dimensions')
    end
end