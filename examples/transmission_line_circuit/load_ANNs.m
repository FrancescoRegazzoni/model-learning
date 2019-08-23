opt.show_history = 0;
opt.idxNet = inf;

ANNmod_3_8 = read_ANNmod_fromfile(problem,'ANNmod_3_8',opt); 
opt.idxNet = 200;
ANNmod_2_8 = read_ANNmod_fromfile(problem,'ANNmod_2_8',opt); 
opt.idxNet = inf;
ANNmod_1_8 = read_ANNmod_fromfile(problem,'ANNmod_1_8',opt); 

ANNmod_3_8_g = read_ANNmod_fromfile(problem,'ANNmod_3_8_g',opt); 
ANNmod_2_8_g = read_ANNmod_fromfile(problem,'ANNmod_2_8_g',opt); 
ANNmod_1_8_g = read_ANNmod_fromfile(problem,'ANNmod_1_8_g',opt);

ANNs_ext = {ANNmod_1_8_g, ANNmod_2_8_g, ANNmod_3_8_g};
ANNs_int = {ANNmod_1_8, ANNmod_2_8, ANNmod_3_8};