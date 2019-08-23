clear

nSVD = 6;
do_save = 0;
    
figure('units','pixel','outerposition',[500 500 700 400]);

for i_fig = [2 1]
    if i_fig==1
        problem = problem_get('thermalblock','TB9_sourceterm.ini');
        sourceterm_pb = 1;
    else
        problem = problem_get('thermalblock','TB9.ini');
        sourceterm_pb = 0;
    end
    problem.goto_dir()
    %HFmod = problem.get_model(problem);

    
    %% ANN model loading
    load_ANNs;

    dataset_def.problem = problem;
    dataset_def.type = 'file';
    dataset_def.source = 'samples_rnd.mat;1:50|samples_x_const_lhs20.mat;:|samples_T10.mat;1:10';
    ds_test = dataset_get(dataset_def);

    models = {ANNmod_g_1_12, ANNmod_g_2_12, ANNmod_g_3_12};
    nmod = length(models);
    for i = 1:nmod
        out= model_compute_error(models{i},ds_test);
        err(i) = out.err_dataset_L2_norm;
    end

    POD_dataset = 'samples_x_const_lhs10.mat;1:10|samples_x_rnd_A.mat;1:50';
    X = build_snapshots_matrix(problem,POD_dataset);

    optPOD.get_full_V = 1;
    [V,outSVD] = POD_projection(X,optPOD);
    %%
    for iN = 1:nSVD
        ep(iN) = POD_getepsilon(outSVD.s,iN);
    end
    %%
    
    if i_fig == 1
        col = 'r';
    else
        col = 'b';
    end

    semilogy(1:nmod, err,[col 'o--'])
    hold on
    semilogy(1:nSVD, ep,[col '+-'])
    pause(1e-16)
end
%%

legend('input conductivity: test error','input conductivity: \epsilon_n', 'input source term: test error','input source term: \epsilon_n')
xlabel('n')
ylabel('err')
grid on
    
if do_save
    create_directory_if_not_found('fig')
    print('fig/TB_sourceterm_errors','-depsc'); 
end