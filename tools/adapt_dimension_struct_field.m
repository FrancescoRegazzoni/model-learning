function output_struct = adapt_dimension_struct_field(input_struct,field_name,dim)
    % Transforms the field 'field_name' of 'input_struct' into a dim-by-1
    % matrix.
    
    if isfield(input_struct,field_name)
        output_struct = setfield(input_struct,field_name, ...
            adapt_dimension(getfield(input_struct,field_name),dim) ...
            );
    else
        output_struct = input_struct;
    end
end