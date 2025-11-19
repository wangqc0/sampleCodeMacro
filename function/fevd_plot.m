function [] = fevd_plot(fevds, fevd_varname, color)
    arguments
        fevds (:, :, :) double
        fevd_varname (1, :) cell = {}
        color (:, 3) double = []
    end
    n_plot_fevd = size(fevds, 3);
    if (numel(fevd_varname) ~= n_plot_fevd) && ~isempty(fevd_varname)
        fevd_varname = strcat({'Variable '}, num2str((1:n_plot_fevd)'))';
    end
    h_max = size(fevds, 1) - 1;
    n_plot_col = floor(sqrt(n_plot_fevd));
    n_plot_row = ceil(n_plot_fevd / n_plot_col);
    tiledlayout(n_plot_row, n_plot_col);
    for i_plot_fevd = 1:n_plot_fevd
        nexttile
        area(0:h_max, fevds(:, :, i_plot_fevd))
        if ~isempty(color)
            colororder(color)
        end
        ylim([0 1])
        title(fevd_varname{i_plot_fevd})
    end
    sgtitle('FEVD')
    leg = legend(fevd_varname, 'Orientation', 'horizontal');
    leg.Layout.Tile = 'south';
end