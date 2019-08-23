function test_solve = model_get_random_test(model)
    % Gets a randomic test for 'model'.
    
    if ~isfield(model,'dt')
        model.dt = model.problem.T * 1e-2;
    end
    test_solve.tt = 0:model.dt:model.problem.T;
    if model.problem.nU > 0
        optRandomU.dim = model.problem.nU;
        optRandomU.umin = model.problem.u_min;
        optRandomU.umax = model.problem.u_max;
        test_solve.uu = get_random_time_course(test_solve.tt,optRandomU);
    end
    
end