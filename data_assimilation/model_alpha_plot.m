function alphafig = model_alpha_plot(mod_learned)

    alpha_original_available = isfield(mod_learned, 'alpha_original');
    dataset_train = getfield(datasetcouple_get(mod_learned.datasets_def), 'train');
    
    if mod_learned.nA == 1
        alphafig = figure();        
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
            xval = mod_learned.alpha_learned(1,iS);
            if alpha_original_available
                yval = mod_learned.alpha_original(1,iS);                    
            else
                yval = xval;
            end
            hhh = plot(xval, yval, 'o');
            set(hhh, 'MarkerFaceColor', get(hhh,'Color')); 
            hold on
            if isfield(dataset_train{iS}, 'label')
                label = dataset_train{iS}.label;
            else
                label = sprintf('test %d', iS);
            end
            text(xval, yval, label)
        end
        if alpha_original_available
            xlabel('\alpha learned')
            ylabel('\alpha original')
        end
    elseif mod_learned.nA == 2 || mod_learned.nA == 3
        alphafig = figure();        
        for i = 1:size(mod_learned.alpha_learned,2)
            if isfield(dataset_train{i}, 'label')
                label = dataset_train{i}.label;
            else
                label = sprintf('test %d', i);
            end
            
            if alpha_original_available
                subplot(1,2,1)
            end
            if mod_learned.nA == 2
                hhh = plot(mod_learned.alpha_learned(1,i),mod_learned.alpha_learned(2,i),'o'); hold on
                text(mod_learned.alpha_learned(1,i),mod_learned.alpha_learned(2,i),label)
            elseif mod_learned.nA == 3
                hhh = plot3(mod_learned.alpha_learned(1,i),mod_learned.alpha_learned(2,i),mod_learned.alpha_learned(3,i),'o'); hold on
                text(mod_learned.alpha_learned(1,i),mod_learned.alpha_learned(2,i),mod_learned.alpha_learned(3,i),label)
            end
            set(hhh, 'MarkerFaceColor', get(hhh,'Color')); 
            
            if alpha_original_available
                title('learned')
            end
            
            if alpha_original_available
                subplot(1,2,2)
                if mod_learned.nA == 2
                    hhh = plot(mod_learned.alpha_original(1,i),mod_learned.alpha_original(2,i),'o'); hold on
                    text(mod_learned.alpha_original(1,i),mod_learned.alpha_original(2,i),label)
                elseif mod_learned.nA == 3
                    hhh = plot(mod_learned.alpha_original(1,i),mod_learned.alpha_original(2,i),mod_learned.alpha_original(3,i),'o'); hold on
                    text(mod_learned.alpha_original(1,i),mod_learned.alpha_original(2,i),mod_learned.alpha_original(3,i),label)
                end
                set(hhh, 'MarkerFaceColor', get(hhh,'Color')); 
                title('original')
            end
        end
    else
        warning('alpha plot is available only if nA = 1')
    end

end