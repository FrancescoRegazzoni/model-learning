function vv_out = interp_time_series(tt_in,vv_in,tt_out,opt)
    % Interpolates the time series (tt_in, vv_in) on the times tt_out.

    opt.dummy = 0;
    if ~isfield(opt,'mode')
        opt.mode = 'pointwise';
    end
    
    is_equal = 0;
    if size(tt_in,2) == size(tt_out,2)
        if norm(tt_in-tt_out)/norm(tt_in) < 1e-14
            is_equal = 1;
        end
    end
    
    if is_equal
        vv_out = vv_in;
    else
        vv_out = zeros(size(vv_in,1),size(tt_out,2));
        switch opt.mode
            case 'pointwise'
                approach = 1;
            case 'mean_forward'
                if length(tt_in) < length(tt_out)
                    % interpolation of a finer grid
                    approach = 1;
                else
                    % interpolation on a coarser grid
                    approach = 2;
                end
            otherwise
                error('interpolation model %s non existent.',opt.mode)
        end
                
        switch approach
            case 1
                for i = 1:size(vv_in,1)
                    vv_out(i,:) = interp1(tt_in,vv_in(i,:),tt_out);
                end
            case 2
                tt_out_ext = [tt_out 2*tt_out(end)-tt_out(end-1)];
                for iT = 1:size(tt_out,2)
                    iT_start = find(tt_in>=tt_out(iT),1,'first');
                    iT_end = find(tt_in<tt_out_ext(iT+1),1,'last');
                    vv_out(:,iT) = mean(vv_in(:,iT_start:iT_end),2);
                end
        end
    end
    
end