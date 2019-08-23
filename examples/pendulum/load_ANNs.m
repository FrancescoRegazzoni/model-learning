opt.show_history = 0;
opt.get_history = 1;

%% load models for phase space test
iANNmod = 1;
ANNmodels_N2{iANNmod} = read_ANNmod_fromfile(problem,'ANNmodels_N2_a',opt); iANNmod = iANNmod +1;
ANNmodels_N2{iANNmod} = read_ANNmod_fromfile(problem,'ANNmodels_N2_b',opt); iANNmod = iANNmod +1;
ANNmodels_N2{iANNmod} = read_ANNmod_fromfile(problem,'ANNmodels_N2_c',opt); iANNmod = iANNmod +1;
ANNmodels_N2{iANNmod} = read_ANNmod_fromfile(problem,'ANNmodels_N2_d',opt); iANNmod = iANNmod +1;

%% load best models
ANNmodels_best{1} = read_ANNmod_fromfile(problem,'ANNmodels_best_N1',opt);
ANNmodels_best{2} = read_ANNmod_fromfile(problem,'ANNmodels_best_N2',opt);