function [dy_dv_bs, dy_de_bs, P_bs] = var_irf_chol_bs(y, Pi, e, h_max, T_burnin, n_bs, options)
    % The function involves a random draw from a uniform distribution and a
    % datasample. Set seed for replicability.
    arguments
        y (:, :) {mustBeNumeric}
        Pi (:, :) {mustBeNumeric}
        e (:, :) {mustBeNumeric}
        h_max (1, 1) {mustBeInteger}
        T_burnin (1, 1) {mustBeInteger}
        n_bs (1, 1) {mustBeInteger}
        options.drift logical = false
    end
    % derived parameters:
    var_e = cov(e');
    P = chol(var_e)';
    n = size(var_e, 1);
    if options.drift
        p = (size(Pi, 2) - 2) / n;
    else
        p = (size(Pi, 2) - 1) / n;
    end
    T = size(y, 2);
    % reshuffle residual
    tau = round((p - .5) + rand(n_bs, 1) * ((T + .5) - (p - .5)));
    % SLIGHTLY FASTER:
    y_initial = cell2mat(arrayfun(@(tau_b) y(:, (tau_b - p + 1):tau_b), tau, 'UniformOutput', false));
    y_initial = pagetranspose(reshape(y_initial', p, n, n_bs));
    e_reshuffle = reshape(e(:, datasample(1:(T - p), (T + T_burnin) * n_bs)), n, T + T_burnin, n_bs);
%     y_initial = nan(n, p, n_bs);
%     e_reshuffle = nan(size(e, 1), size(e, 2) + T_burnin + p, n_bs);
%     for i_bs = 1:n_bs
%         y_initial(:, :, i_bs) = y(:, (tau(i_bs) - p + 1):tau(i_bs));
%         e_reshuffle(:, :, i_bs) = e(:, datasample(1:(T - p), T + T_burnin));
%     end
    % generate derived y series
    % FASTER:
    y_new = nan(n, T + T_burnin, n_bs);
    y_new(:, 1:p, :) = y_initial;
    for t = (p + 1):(T + T_burnin)
        if options.drift
            x = nan(2 + n * p, 1, n_bs);
            x(1, 1, :) = 1;
            x(2, 1, :) = t - p - 1 - T_burnin;
            for p_i = 1:p
                x((2 + ((n * (p_i - 1) + 1):(n * p_i))), 1, :) = y_new(:, (t - p_i), :);
            end
        else
            x = nan(1 + n * p, 1, n_bs);
            x(1, 1, :) = 1;
            for p_i = 1:p
                x((1 + ((n * (p_i - 1) + 1):(n * p_i))), 1, :) = y_new(:, (t - p_i), :);
            end
        end
        y_new(:, t, :) = pagemtimes(Pi, x) + e_reshuffle(:, t, :);
    end
%     y_new = nan(size(y, 1), size(y, 2) + T_burnin, n_bs);
%     for i_bs = 1:n_bs
%         y_new_i = y_new(:, :, i_bs);
%         y_new_i(:, 1:p) = y_initial(:, :, i_bs);
%         for t = (p + 1):(T + T_burnin)
%             x = zeros(1 + n * p, 1);
%             x(1) = 1;
%             for p_i = 1:p
%                 x((1 + ((n * (p_i - 1) + 1):(n * p_i)))) = y_new_i(:, (t - p_i));
%             end
%             y_new_i(:, t) = Pi * x + e_reshuffle(:, t, i_bs);
%         end
%         y_new(:, :, i_bs) = y_new_i;
%     end
    y_new = y_new(:, (T_burnin + 1):end, :);
    % obtain Pi and Omega
    Pi_bs = nan(size(Pi, 1), size(Pi, 2), n_bs);
    Omega_bs = nan(size(var_e, 1), size(var_e, 2), n_bs);
    P_bs = nan(size(P, 1), size(P, 2), n_bs);
    for i_bs = 1:n_bs
        if options.drift
            [Pi_bs_i, e_i] = var_ols(y_new(:, :, i_bs), p, 'drift', true);
        else
            [Pi_bs_i, e_i] = var_ols(y_new(:, :, i_bs), p);
        end
        var_e_i = cov(e_i');
        P_bs_i = chol(var_e_i)';
        Pi_bs(:, :, i_bs) = Pi_bs_i;
        Omega_bs(:, :, i_bs) = var_e_i;
        P_bs(:, :, i_bs) = P_bs_i;
    end
    % obtain Psi
    % SLIGHTLY FASTER
    F_bs = zeros(n * p, n * p, n_bs);
    if options.drift
        F_bs(1:n, :, :) = Pi_bs(:, 3:end, :);
    else
        F_bs(1:n, :, :) = Pi_bs(:, 2:end, :);
    end
    F_bs((n + 1):end, 1:(end - n), :) = repmat(eye(n * (p - 1), n * (p - 1)), 1, 1, n_bs);
    %F_bs_h = reshape(repmat(F_bs, 1, 1, h_max), n * p, n * p, n_bs, h_max);
    %F_bs_h = permute(F_bs_h, [1 2 4 3]);
%     F_bs = zeros(n * p, n * p, n_bs);
%     F_bs_h = nan(n * p, n * p, h_max, n_bs);
%     F_bs(1:n, :, :) = Pi_bs(:, 2:end, :);
%     for i_bs = 1:n_bs
%         F_bs((n + 1):end, 1:(end - n), i_bs) = eye(n * (p - 1), n * (p - 1));
%         F_bs_h(:, :, :, i_bs) = repmat(F_bs(:, :, i_bs), 1, 1, h_max);
%     end
    % FASTER
    F_bs_h = nan(n * p, n * p, h_max, n_bs);
    for i_bs = 1:n_bs
        F_i = F_bs(:, :, i_bs);
        F_h_i = nan(n * p, n * p, h_max);
        F_h_i_h = eye(n * p);
        for h = 1:h_max
            F_h_i_h = F_h_i_h * F_i;
            F_h_i(:, :, h) = F_h_i_h;
        end
        F_bs_h(:, :, :, i_bs) = F_h_i;
    end
    Psi_bs = F_bs_h(1:n, 1:n, :, :);
    Psi_bs = cat(3, repmat(eye(n), 1, 1, 1, n_bs), Psi_bs);
%     Psi_bs = nan(n, n, h_max + 1, n_bs);
%     for i_bs = 1:n_bs
%         F_h_i = F_bs_h(:, :, :, i_bs);
%         Psi_i = nan(n, n, h_max);
%         for h = 1:h_max
%             F_h_i(:, :, h) = F_h_i(:, :, h) ^ h;
%             Psi_i(:, :, h) = F_h_i(1:n, 1:n, h);
%         end
%         % add t = 0
%         Psi_i = cat(3, eye(n), Psi_i);
%         F_bs_h(:, :, :, i_bs) = F_h_i;
%         Psi_bs(:, :, :, i_bs) = Psi_i;
%     end
    dy_de_bs = Psi_bs;
    % muitiply Psi and P
    % FASTER
    Psi_s_bs = pagemtimes(Psi_bs, permute(reshape(repmat(P_bs, 1, 1, h_max + 1), n, n, n_bs, h_max + 1), [1 2 4 3]));
%     Psi_s_bs = nan(size(Psi_bs));
%     for i_bs = 1:n_bs
%         Psi_i = Psi_bs(:, :, :, i_bs);
%         P_i = P_bs(:, :, i_bs);
%         Psi_s_i = pagemtimes(Psi_i, P_i);
%         Psi_s_bs(:, :, :, i_bs) = Psi_s_i;
%     end
    dy_dv_bs = Psi_s_bs;
end