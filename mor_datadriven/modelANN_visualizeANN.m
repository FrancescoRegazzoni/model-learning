function modelANN_visualizeANN(model)

figure()

opts.plot_labels_IO = 1;
for i = 1:model.problem.nU
    opts.label_in{i} = sprintf('{u}_{%d}',i);
end
for i = 1:model.numnF(1)-model.problem.nU
    iPoint = model.problem.nU+i;
    opts.label_in{iPoint} = sprintf('{x}_{%d}',i);
end

for i = 1:model.numnF(end)
    opts.label_out{i} = sprintf('\\dot{x}_{%d}',i);
end

ANN_visualize(model.numnF,model.wF,model.thetaF,opts)

end