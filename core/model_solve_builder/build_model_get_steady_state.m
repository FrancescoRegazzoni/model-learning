function build_model_get_steady_state()

    advance_types = get_advance_types();
    
    %% opening file
    file_out = fopen('model_get_steady_state.m','w');
    
    %% writing header
    write_copy_file(file_out,'model_get_steady_state_header.m')            
    fprintf(file_out,'\tswitch model.advance_type\n');    
    for i_adv = 1:length(advance_types)
        fprintf(file_out,'\t\tcase ''%s''\n',advance_types{i_adv});        
        write_copy_file(file_out,'model_get_steady_state_loop_header.m')
        advance_types_print_switch(file_out,'\t\t\t',advance_types{i_adv}, ...
                    struct('u_name','u_steady','dt_name','dt'));
        write_copy_file(file_out,'model_get_steady_state_loop_footer.m')         
    end      
    fprintf(file_out,'\t\totherwise\n');
    fprintf(file_out,'\t\t\terror(''unknown advance type %%s'',model.advance_type)\n');
    fprintf(file_out,'\tend\n');
    fprintf(file_out,'end\n');
    
    %% closing file
    fclose(file_out);
    
    disp("Script generated! Do not forget to copy model_get_steady_state.m into the 'core' folder if you want to keep changes.")
end