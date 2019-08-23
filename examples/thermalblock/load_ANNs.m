opt.show_history = 0;
opt.idxNet = inf;

if sourceterm_pb
    ANNmod_3_12 = read_ANNmod_fromfile(problem,'ANNmod_3_12',opt);
    ANNs_int = {ANNmod_3_12};  
    
    ANNmod_g_3_12 = read_ANNmod_fromfile(problem,'ANNmod_g_3_12',opt); 
    ANNmod_g_2_12 = read_ANNmod_fromfile(problem,'ANNmod_g_2_12',opt); 
    ANNmod_g_1_12 = read_ANNmod_fromfile(problem,'ANNmod_g_1_12',opt);
    ANNs_ext = {ANNmod_g_1_12, ANNmod_g_2_12, ANNmod_g_3_12};
else
    ANNmod_5_16 = read_ANNmod_fromfile(problem,'ANNmod_5_16',opt); 
    ANNmod_4_12 = read_ANNmod_fromfile(problem,'ANNmod_4_12',opt); 
    ANNmod_3_12 = read_ANNmod_fromfile(problem,'ANNmod_3_12',opt); 
    ANNs_int = {ANNmod_3_12, ANNmod_4_12, ANNmod_5_16};    
    
    ANNmod_g_4_12 = read_ANNmod_fromfile(problem,'ANNmod_g_4_12',opt); 
    ANNmod_g_3_12 = read_ANNmod_fromfile(problem,'ANNmod_g_3_12',opt); 
    ANNmod_g_2_12 = read_ANNmod_fromfile(problem,'ANNmod_g_2_12',opt); 
    ANNmod_g_1_12 = read_ANNmod_fromfile(problem,'ANNmod_g_1_12',opt);    
    ANNs_ext = {ANNmod_g_1_12, ANNmod_g_2_12, ANNmod_g_3_12, ANNmod_g_4_12};
end
