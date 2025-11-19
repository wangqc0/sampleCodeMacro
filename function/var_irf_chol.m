function [dy_dv, dy_de, P] = var_irf_chol(Pi, e, h_max, options)
    arguments
        Pi (:, :) {mustBeNumeric}
        e (:, :) {mustBeNumeric}
        h_max (1, 1) {mustBeInteger}
        options.drift logical = false
    end
    var_e = cov(e');
    P = chol(var_e)';
    n = size(var_e, 1);
    if options.drift
        p = (size(Pi, 2) - 2) / n;
    else
        p = (size(Pi, 2) - 1) / n;
    end
    % use the companion form to get Psi from Phi
    F = zeros(n * p, n * p);
    if options.drift
        F(1:n, :) = Pi(:, 3:end);
    else
        F(1:n, :) = Pi(:, 2:end);
    end
    F((n + 1):end, 1:(end - n)) = eye(n * (p - 1));
    % NOT SUBSTANTIALLY FASTER
%     F_raise_power = @(x, y) (x ^ y);
%     F_extract_upperleft = @(x, y) (x(1:y, 1:y));
%     Psi = cell2mat(arrayfun(@(h) F_extract_upperleft(F_raise_power(F, h), n), 1:h_max, 'UniformOutput', false));
%     Psi = reshape(Psi, n, n, h_max);
    F_h = repmat(F, 1, 1, h_max);
    Psi = nan(n, n, h_max);
    for h = 1:h_max
        F_h(:, :, h) = F_h(:, :, h) ^ h;
        Psi(:, :, h) = F_h(1:n, 1:n, h);
    end
    % add t = 0
    Psi = cat(3, eye(n), Psi);
    dy_de = Psi;
    dy_dv = pagemtimes(dy_de, P);
end
