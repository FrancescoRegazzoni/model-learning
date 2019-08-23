opt.show_history = 0;
opt.get_history = 1;

ANNmod_2 = read_model_fromfile(problem,'ANN_N2',opt); 
ANNmod_4 = read_model_fromfile(problem,'ANN_N4',opt); 
ANNmod_6 = read_model_fromfile(problem,'ANN_N6',opt);
ANNmod_8 = read_model_fromfile(problem,'ANN_N8',opt);

ANNs = {ANNmod_2,ANNmod_4,ANNmod_6,ANNmod_8};