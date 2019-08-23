function dataset_plot_legend(lines,filename)

    cmap = get(0, 'DefaultAxesColorOrder');
    % col_exact       = cmap(1,:);
    % col_y           = cmap(2,:);
    % col_u1          = cmap(3,:);
    % col_u2          = cmap(4,:);
    % sty_exact       = '-';
    % sty_y           = '--';
    % sty_u1          = ':';
    % sty_u2          = ':';

    figure('units','pixel','outerposition',[500 500 100 200]);
    for i=1:length(lines)
        if ~isfield(lines{i},'style')
            lines{i}.style = '-';
        end
        if ~isfield(lines{i},'color')
            lines{i}.color = cmap(i,:);
        end
        if ~isfield(lines{i},'width')
            lines{i}.width = 1;
        end
        plot(0,0,lines{i}.style,'linewidth',lines{i}.width,'Color',lines{i}.color); 
        hold on;
        names{i} = lines{i}.name;
    end
    legend(names)  
    set(gca, 'visible', 'off');

    print(filename,'-depsc'); 
end