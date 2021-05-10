function fig_alpha_plot(mod_learned)
    
    if isfield(mod_learned,'alpha_to_alpha')
        npt = 1e2;
        aa = linspace(min(mod_learned.alpha_original),max(mod_learned.alpha_original),npt);
        for ia = 1:npt
            aa_tilde(ia) = mod_learned.alpha_to_alpha(aa(ia));
        end
        plot(aa_tilde,aa,'k-')
        hold on
    end
    for iS = 1:length(mod_learned.alpha_learned)
        hhh = plot(mod_learned.alpha_learned(1,iS), mod_learned.alpha_original(1,iS),'o');
        set(hhh, 'MarkerFaceColor', get(hhh,'Color')); 
        hold on
    end
    if isfield(mod_learned,'alpha_to_alpha')
        axis([0 1.1 0 1])
%         daspect([1.8 1 1])
    else
        axis([-1 1 0 1])
        daspect([1.8 1 1])
    end
    xlabel('$\hat\alpha$ (learned)', 'interpreter', 'latex')
    ylabel('$\alpha$ (original)', 'interpreter', 'latex')
%     set(gcf,'position',figure_size_alpha)
%     aa = get(gcf,'position')
end