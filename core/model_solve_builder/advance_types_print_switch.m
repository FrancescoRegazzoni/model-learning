function advance_types_print_switch(file_out,prefix,adv_type,opt)
    opt.dummy = 0;
    if ~isfield(opt,'u_name')
        opt.u_name = 'u';
    end
    if ~isfield(opt,'dt_name')
        opt.dt_name = 'dt';
    end
    
    file_adv = fopen('model_solve_advance_types.m','r');
    found_adv = 0;
    while ~feof(file_adv)
        line = fgetl(file_adv);
        if startsWith(line,'%%CASE')
            if startsWith(line,sprintf('%%%%CASE ''%s''',adv_type))
                found_adv = 1;
            else
                if found_adv
                    break
                end
            end
        else
            if found_adv
                line = strrep(line,'uu_eff(:,iT)',opt.u_name);
                line = strrep(line,'dtt(iT)',opt.dt_name);
                fprintf(file_out,[prefix '%s\n'],line);
            end
        end
    end
    fclose(file_adv);
end