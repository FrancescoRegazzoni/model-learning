function out = plotting_get_axes(num_rows,num_cols,opt)

    opt.dummy = 0;
    %%%%%%%% size
    if ~isfield(opt,'subplot_width') 
        opt.subplot_width = 120;
    end
    if ~isfield(opt,'subplot_heigth') 
        opt.subplot_heigth = 80;
    end
    %%%%%%%% margins
    if ~isfield(opt,'left_margin') 
        opt.left_margin = 20;
    end
    if ~isfield(opt,'right_margin') 
        opt.right_margin = 20;
    end
    if ~isfield(opt,'bottom_margin') 
        opt.bottom_margin = 20;
    end
    if ~isfield(opt,'top_margin') 
        opt.top_margin = 2;
    end
    if ~isfield(opt,'margin_plot_h') 
        opt.margin_plot_h = 20;
    end
    if ~isfield(opt,'margin_plot_v') 
        opt.margin_plot_v = 20;
    end
    %%%%%%%% legend
    if ~isfield(opt,'plot_legend')
        opt.plot_legend = 0;
    end
    if ~isfield(opt,'legend_orientation')
        opt.legend_orientation = 'h';
    end
    if ~isfield(opt,'legend_width')
        opt.legend_width = 150;
    end
    if ~isfield(opt,'legend_heigth')
        opt.legend_heigth = 20;
    end
    if ~isfield(opt,'legend_margin') 
        opt.legend_margin = 2;
    end
    
    
    out.width = opt.left_margin + opt.right_margin ...
              + num_cols*opt.subplot_width ...
              + (num_cols-1)*opt.margin_plot_h;
    out.heigth = opt.bottom_margin + opt.top_margin ...
              + num_rows*opt.subplot_heigth ...
              + (num_rows-1)*opt.margin_plot_v;
    if opt.plot_legend
        if opt.legend_orientation == 'h'
            out.heigth = out.heigth + opt.legend_heigth + opt.legend_margin;
        else
            out.width = out.width + opt.legend_width + opt.legend_margin;
        end
    end
    out.normalization = [out.width out.heigth out.width out.heigth];
    if opt.plot_legend
        if opt.legend_orientation == 'h'
            out.legend_coord = [(out.width - opt.legend_width)/2, ...
                                opt.legend_margin, ...
                                opt.legend_width, ...
                                opt.legend_heigth];
        else
            out.legend_coord = [out.width - opt.legend_width - opt.legend_margin, ...
                                (out.heigth - opt.legend_heigth)/2, ...
                                opt.legend_width, ...
                                opt.legend_heigth];
            
        end
        out.legend_coord_norm = out.legend_coord ./ out.normalization;
        out.get_legend_coord = @get_legend_coord;
        out.get_legend_coord_norm = @get_legend_coord_norm;
    end
    
    
    row_Curr = 0;
    col_Curr = 0;
    for j = 1:num_cols*num_rows
        leftCurr = opt.left_margin ...
                 + col_Curr * (opt.subplot_width + opt.margin_plot_h);
        bottomCurr = out.heigth ...
                   - opt.top_margin ...
                   - (row_Curr+1)*opt.subplot_heigth ...
                   - row_Curr*opt.margin_plot_v;
        out.axs{j}.coord = [leftCurr , bottomCurr , opt.subplot_width , opt.subplot_heigth];
        out.axs{j}.coord_norm = out.axs{j}.coord ./ out.normalization;
        out.axs{j}.col = col_Curr + 1;
        out.axs{j}.row = row_Curr + 1;
        col_Curr = col_Curr+1;
        if col_Curr == num_cols
            col_Curr = 0;
            row_Curr = row_Curr+1;
        end
    end    
    
    function new_coord = get_legend_coord(coord)
        coord = coord .* out.normalization;
        if opt.legend_orientation == 'h'
            new_coord = [(out.width - coord(3))/2, ...
                     opt.legend_margin, ...
                     coord(3), ...
                     opt.legend_heigth]; 
        else
            new_coord = [out.width - opt.legend_width - opt.legend_margin, ...
                     (out.heigth - coord(4))/2, ...
                     opt.legend_width, ...
                     coord(4)];            
        end
        
    end

    function new_coord = get_legend_coord_norm(coord)
        new_coord = get_legend_coord(coord) ./ out.normalization;
    end
end