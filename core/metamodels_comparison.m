function metamodels_comparison(mod_HF,mod_learned)
    % Compares a learned meta-model with the HF model.

    if mod_HF.problem.nY ~= 1 || mod_HF.nX ~= 1 || mod_HF.nA ~= 1 || mod_learned.nX ~= 1 || mod_learned.nA ~= 1
        warning('wrong dimensions in metamodels')
        return
    end
    if ~isfield(mod_learned,'alpha_to_alpha')
        warning('alpha_to_alpha function not avalable')
        return
    end
    
    nptx = 40;
    npta = 40;
    x_min = mod_HF.problem.y_min;
    x_max = mod_HF.problem.y_max;
    a_min = mod_HF.alpha_min;
    a_max = mod_HF.alpha_max;
    xx = linspace(x_min,x_max,nptx);
    aa = linspace(a_min,a_max,npta);
    [XX,AA] = meshgrid(xx,aa);
    FF_HF = zeros(npta,nptx);
    FF_NN = zeros(npta,nptx);
    for ix = 1:nptx
        for ia = 1:npta
            FF_HF(ia,ix) = mod_HF.f_alpha(xx(ix),[],aa(ia));
            FF_NN(ia,ix) = mod_learned.f_alpha(xx(ix),[],mod_learned.alpha_to_alpha(aa(ia)));
        end
    end

    Fmin = min(min(FF_HF(:)),min(FF_NN(:)));
    Fmax = max(max(FF_HF(:)),max(FF_NN(:)));

    figure('units','normalized','outerposition',[0 0 1 1]);
    ncols = 3;
    nrows = 2;

    subplot(nrows,ncols,1)
    surf(XX,AA,FF_HF)
    title('HF model')
    axis([x_min x_max a_min a_max Fmin Fmax])
    xlabel('x'); ylabel('\alpha'); zlabel('f')

    subplot(nrows,ncols,2)
    surf(XX,AA,FF_NN)
    title('Learned model')
    axis([x_min x_max a_min a_max Fmin Fmax])
    xlabel('x'); ylabel('\alpha'); zlabel('f')

    subplot(nrows,ncols,3)
    surf(FF_HF, 'FaceAlpha', 0.8, 'EdgeColor', 'b', 'LineWidth', 2.0);
    hold on
    surf(FF_NN, 'FaceAlpha', 0.2, 'EdgeColor', 'r', 'LineWidth', 0.2);
    xlabel('x'); ylabel('\alpha'); zlabel('f')
    hold off

    subplot(nrows,ncols,4)
    surf(XX,AA,abs(FF_HF-FF_NN))
    title('error')
    xlabel('x'); ylabel('\alpha');

    subplot(nrows,ncols,5)
    surf(XX,AA,log10(abs(FF_HF-FF_NN)))
    title('log10 error')
    xlabel('x'); ylabel('\alpha');

end