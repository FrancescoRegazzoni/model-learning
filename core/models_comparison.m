function models_comparison(mod1,mod2)
    % Compares two models (if the dimensions allow to do a comparison plot).

    if (mod1.problem.samples_variability && ~mod1.problem.particularized) || ...
       (mod2.problem.samples_variability && ~mod2.problem.particularized) 
        warning('models with samples variability should be compared thorugh the function metamodels_comparison')
        return
    end
    
    if mod1.problem.nU ~= 1 || mod1.nX ~= 1 || mod2.problem.nU ~= 1 || mod2.nX ~= 1
        warning('wrong dimensions in models')
        return
    end
    
    nptx = 40;
    nptu = 40;
    u_min = min(mod1.problem.u_min, mod2.problem.u_min);
    u_max = max(mod1.problem.u_max, mod2.problem.u_max);
    if isequal(mod1.output_type, 'insidestate') && isequal(mod2.output_type, 'insidestate')
        x_min = min(mod1.problem.y_min, mod2.problem.y_min);
        x_max = max(mod1.problem.y_max, mod2.problem.y_max);
    else
        error('not (yet) implemented')
    end
    xx = linspace(x_min,x_max,nptx);
    uu = linspace(u_min,u_max,nptu);
    [XX,UU] = meshgrid(xx,uu);
    FF_1 = zeros(nptu,nptx);
    FF_2 = zeros(nptu,nptx);
    for ix = 1:nptx
        for iu = 1:nptu
            FF_1(iu,ix) = mod1.f(xx(ix),uu(iu));
            FF_2(iu,ix) = mod2.f(xx(ix),uu(iu));
        end
    end

    Fmin = min(min(FF_1(:)),min(FF_2(:)));
    Fmax = max(max(FF_1(:)),max(FF_2(:)));
    
    norm_F = (Fmax-Fmin);
    relErr = abs(FF_1-FF_2)/norm_F;

    ncols = 3;
    nrows = 2;

    subplot(nrows,ncols,1)
    surf(XX,UU,FF_1)
    title('HF model')
    axis([x_min x_max u_min u_max Fmin Fmax])
    xlabel('x'); ylabel('u'); zlabel('f')

    subplot(nrows,ncols,2)
    surf(XX,UU,FF_2)
    title('Learned model')
    axis([x_min x_max u_min u_max Fmin Fmax])
    xlabel('x'); ylabel('u'); zlabel('f')

    subplot(nrows,ncols,3)
    surf(FF_1, 'FaceAlpha', 0.8, 'EdgeColor', 'b', 'LineWidth', 2.0);
    hold on
    surf(FF_2, 'FaceAlpha', 0.2, 'EdgeColor', 'r', 'LineWidth', 0.2);
    xlabel('x'); ylabel('u'); zlabel('f')
    hold off

    subplot(nrows,ncols,4)
    surf(XX,UU,relErr)
    title('relative error')
    xlabel('x'); ylabel('u');

    subplot(nrows,ncols,5)
    surf(XX,UU,log10(relErr))
    title('log10 relative error')
    xlabel('x'); ylabel('u');

end