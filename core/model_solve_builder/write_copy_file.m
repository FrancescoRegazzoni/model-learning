function write_copy_file(file_out,name_file_in)
    
    file_in = fopen(name_file_in,'r');
    
    %% writing header
    while ~feof(file_in)
        fprintf(file_out,'%s\n',fgetl(file_in));
    end
    
    fclose(file_in);
    
end