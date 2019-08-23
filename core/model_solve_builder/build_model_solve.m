function build_model_solve()

    advance_types = get_advance_types();
    
    %% opening file
    file_out = fopen('model_solve.m','w');
    
    %% writing header
    write_copy_file(file_out,'model_solve_header.m')
            
    %% time loop
    fprintf(file_out,'\t\tswitch model.output_type\n');
    for i_output_type = 1:3
        soutput = 'yy(:,iT - idx_tt_y_first + 1) = ';
        switch i_output_type
            case 1
                fprintf(file_out,'\t\tcase ''insidestate''\n');
                soutput = [soutput 'x(1:nY);'];
            case 2
                fprintf(file_out,'\t\tcase ''nonlinear''\n');
                soutput = [soutput 'model.g(x);'];
            case 3
                fprintf(file_out,'\t\tcase ''linear''\n');
                soutput = [soutput 'model.G*x+model.g0;'];                
        end        
        for i_save_x = 1:3
            switch i_save_x
                case 1
                    fprintf(file_out,'\t\t\tif opt.do_plot_x || (opt.save_x && opt.save_x_freq == 1 )\n');
                    s_output_append = {'xx  = [xx x];'};
                case 2
                    fprintf(file_out,'\t\t\telseif opt.save_x\n');
                    s_output_append = {'if mod(iT,opt.save_x_freq) == 0','    xx  = [xx x];','end'};
                case 3
                    fprintf(file_out,'\t\t\telse\n');
                    s_output_append = {};
            end                    
            fprintf(file_out,'\t\t\t\tswitch model.advance_type\n');
            for i_adv = 1:length(advance_types)
                fprintf(file_out,'\t\t\t\t\tcase ''%s''\n',advance_types{i_adv});
                fprintf(file_out,'\t\t\t\t\t\tfor iT = 2:nT\n');
                %model_solve_advance_types(advance_types{i_adv},file_out)
                advance_types_print_switch(file_out,'\t\t\t\t\t\t',advance_types{i_adv}, ...
                    struct('u_name','uu_eff(:,iT)','dt_name','dtt(iT)'));
                % compute output y
                fprintf(file_out,'\t\t\t\t\t\t\t%s\n','if iT >= idx_tt_y_first && iT<= idx_tt_y_last');
                fprintf(file_out,'\t\t\t\t\t\t\t\t%s\n',soutput);  
                fprintf(file_out,'\t\t\t\t\t\t\t%s\n','end');
               
                % output append
                for i_out_app = 1:length(s_output_append)
                    fprintf(file_out,'\t\t\t\t\t\t\t%s\n',s_output_append{i_out_app});
                end
                
                fprintf(file_out,'\t\t\t\t\t\tend\n');
            end
            fprintf(file_out,'\t\t\t\t\totherwise\n');
            fprintf(file_out,'\t\t\t\t\t\terror(''unknown advance type %%s'',model.advance_type)\n');
            fprintf(file_out,'\t\t\t\tend\n');
        end
        fprintf(file_out,'\t\t\tend\n');  
    end
    fprintf(file_out,'\t\totherwise\n'); 
    fprintf(file_out,'\t\t\t%s\n','error(''unknown output type %s'',model.output_type)'); 
    fprintf(file_out,'\t\tend\n');    
    
    %% writing footer
    write_copy_file(file_out,'model_solve_footer.m')
    
    %% closing file
    fclose(file_out);
    
    disp("Script generated! Do not forget to copy model_solve.m into the 'core' folder if you want to keep changes.")
end