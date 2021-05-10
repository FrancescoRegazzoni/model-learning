function modelclass = modelclass_ANN_alphaconstrained(optionsfile,problem,N,N_alpha,useG)
    modelclass = modelclass_ANN(optionsfile,problem,N,N_alpha,useG,1);
end