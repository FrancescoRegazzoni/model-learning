function y = model_get_output(model,x)
    switch model.output_type
        case 'insidestate'
            y = x(1:model.problem.nY);
        case 'nonlinear'
            y = model.g(x);
        case 'linear'
            y = model.G*x+model.g0;
    end
end