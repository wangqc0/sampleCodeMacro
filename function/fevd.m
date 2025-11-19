function fevds = fevd(Psi, P, ind_cum, Delta)
    % output: page is response variable, column is shock variable, row is
    % time lag
    arguments
        Psi (:, :, :) double % dy_de
        P (:, :) double
        ind_cum (:, 1) {mustBeInteger} = []
        Delta (:, :) double = eye(size(P)) % For VAR results from Cholesky decomposition, identity matrix works
    end
    if (size(Psi, 1) ~= size(Psi, 2)) || (size(P, 1) ~= size(P, 2)) || (size(Delta, 1) ~= size(Delta, 2))
        error('Psi, P and Delta must be squared matrices')
    end
    if (size(Psi, 1) ~= size(P, 1)) || (size(P, 1) ~= size(Delta, 1))
        error('Psi, P and Delta must have the same size')
    end
    n_var = size(Psi, 1);
    n_h = size(Psi, 3);
    var_fe_h = pagemtimes(pagemtimes(Psi, P * Delta), pagemtimes(P', pagetranspose(Psi)));
    var_fe = cumsum(var_fe_h, 3);
    var_fe_shocks = nan(n_var, n_var, n_h, n_var);
    for i_shock = 1:n_var
        Delta_i = zeros(n_var);
        Delta_i(i_shock, i_shock) = Delta(i_shock, i_shock);
        var_fe_h_i = pagemtimes(pagemtimes(Psi, P * Delta_i), pagemtimes(P', pagetranspose(Psi)));
        var_fe_i = cumsum(var_fe_h_i, 3);
        var_fe_shocks(:, :, :, i_shock) = var_fe_i;
    end
    if ~isempty(ind_cum)
        Psi_cum = cumsum(Psi, 3);
        var_fe_h_cum = pagemtimes(pagemtimes(Psi_cum, P * Delta), pagemtimes(P', pagetranspose(Psi_cum)));
        var_fe_cum = cumsum(var_fe_h_cum, 3);
        var_fe_cum_shocks = nan(n_var, n_var, n_h, n_var);
        for i_shock = 1:n_var
            Delta_i = zeros(n_var);
            Delta_i(i_shock, i_shock) = Delta(i_shock, i_shock);
            var_fe_cum_h_i = pagemtimes(pagemtimes(Psi_cum, P * Delta_i), pagemtimes(P', pagetranspose(Psi_cum)));
            var_fe_cum_i = cumsum(var_fe_cum_h_i, 3);
            var_fe_cum_shocks(:, :, :, i_shock) = var_fe_cum_i;
        end
    end
    fevds = nan(n_h, n_var, n_var);
    for i_response = 1:n_var
        if any(i_response == ind_cum)
            fevd_i = reshape(var_fe_cum_shocks(i_response, i_response, :, :), n_h, n_var, 1, 1) ./ reshape(var_fe_cum(i_response, i_response, :), n_h, 1, 1, 1);
        else
            fevd_i = reshape(var_fe_shocks(i_response, i_response, :, :), n_h, n_var, 1, 1) ./ reshape(var_fe(i_response, i_response, :), n_h, 1, 1, 1);
        end
        fevds(:, :, i_response) = fevd_i;
    end
end