function [] = var_irf_plot(dy_dv, ind_shock, size_shock, ind_cum, dy_dv_bs, bias_correct, irf_varname, irf_band_confidence)
    arguments
        dy_dv (:, :, :) double
        ind_shock (1, 1) {mustBeInteger}
        size_shock (1, 1) double = 1
        ind_cum (:, 1) {mustBeInteger} = []
        dy_dv_bs (:, :, :, :) double = []
        bias_correct (1, 1) logical = 0
        irf_varname (1, :) cell = {}
        irf_band_confidence (1, 1) {mustBeInRange(irf_band_confidence, 60, 100, 'exclude-upper')} = 90
    end
    if size(dy_dv, 1) ~= size(dy_dv, 2)
        error(['In the dy/dv array, the first index is the response variable, ' ...
            'the second index is the shock variable, and they must be equivalent'])
    end
    if size(dy_dv_bs, 1) ~= size(dy_dv_bs, 2)
        error(['In the bootstrapped dy/dv array, the first index is the response variable, ' ...
            'the second index is the shock variable, and they must be equivalent'])
    end
    n_plot_irf = size(dy_dv, 1);
    if (numel(irf_varname) ~= n_plot_irf) && ~isempty(irf_varname)
        irf_varname = strcat({'Variable '}, num2str((1:n_plot_irf)'))';
    end
    h_max = size(dy_dv, 3) - 1;
    irf_response_shock = size_shock * permute(dy_dv(:, ind_shock, :), [3, 1, 2]);
    if ~isempty(ind_cum)
        for index_i_response = 1:numel(ind_cum)
            i_response = ind_cum(index_i_response, 1);
            irf_response_shock(:, i_response) = cumsum(irf_response_shock(:, i_response), 1);
        end
    end
    if ~isempty(dy_dv_bs)
        if (size(dy_dv_bs, 1) ~= size(dy_dv_bs, 1))
            error(['The number of variable in bootstrapped dy/dv should match the number ' ...
                'in dy/dv'])
        end
        if (size(dy_dv_bs, 3) ~= size(dy_dv_bs, 3))
            error(['The time lags of bootstrapped dy/dv should match the time lags ' ...
                'of dy/dv'])
        end
        irf_response_shock_bs = size_shock * permute(dy_dv_bs(:, ind_shock, :, :), [3 1 4 2]);
        if ~isempty(ind_cum)
            for index_i_response = 1:numel(ind_cum)
                i_response = ind_cum(index_i_response, 1);
                irf_response_shock_bs(:, i_response, :) = cumsum(irf_response_shock_bs(:, i_response, :), 1);
            end
        end
        irf_response_shock_bs_upper = prctile(irf_response_shock_bs, 100 - (100 - irf_band_confidence) / 2, 3);
        irf_response_shock_bs_lower = prctile(irf_response_shock_bs, (100 - irf_band_confidence) / 2, 3);
        irf_response_shock_bs_mean = mean(irf_response_shock_bs, 3);
    end
    n_plot_col = floor(sqrt(n_plot_irf));
    n_plot_row = ceil(n_plot_irf / n_plot_col);
    for i_plot_irf = 1:n_plot_irf
        subplot(n_plot_row, n_plot_col, i_plot_irf)
        plot(0:h_max, irf_response_shock(:, i_plot_irf), 'k-')
        if exist('irf_response_shock_bs', 'var')
            hold on
            if bias_correct
                plot(0:h_max, [irf_response_shock_bs_lower(:, i_plot_irf) irf_response_shock_bs_upper(:, i_plot_irf)] + (irf_response_shock(:, i_plot_irf) - irf_response_shock_bs_mean(:, i_plot_irf)), 'k--')
            else
                plot(0:h_max, [irf_response_shock_bs_lower(:, i_plot_irf) irf_response_shock_bs_upper(:, i_plot_irf)], 'k--')
            end
        end
        if ~isempty(irf_varname)
            title(irf_varname{i_plot_irf})
        end
    end
end