function v_out = renormalize_m1p1(v_in,v_min,v_max)
    % Maps 'v_in' from the interval (v_min,v_max) to the interval (-1,1).
    
    v_out = (2*v_in - v_min - v_max)./(v_max-v_min);
end