function v_out = renormalize_back_m1p1(v_in,v_min,v_max)
    % Maps 'v_in' from the interval (-1,1) to the interval (v_min,v_max).
    
    v_out = (v_in.*(v_max-v_min) + v_min + v_max)/2;
end