function model_show_example(model)
    % Performs a randomic test for 'model' and plots the result.

    if model.problem.metaproblem
        model = metamodel_particularize(model);
    end

    test_solve = model_get_random_test(model);
    opt_solve.do_plot = 1;

    figure();
    model_solve(test_solve,model,opt_solve);
end