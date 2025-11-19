function [Pi, e, k, X] = var_ols(y, p, options)
    arguments
        y (:, :) {mustBeNumeric}
        p (1, 1) {mustBeInteger}
        options.drift (1, 1) logical = false
    end
    % y: (n * T), place the earliest value at the beginning
    [n, T] = size(y);
    Y = y(:, (p + 1):T);
    if options.drift
        X = zeros(2 + n * p, T - p);
        X(2, :) = 0:(T - p - 1);
    else
        X = zeros(1 + n * p, T - p);
    end
    X(1, :) = 1;
    for p_i = 1:p
        % y_{t - p_i}
        if options.drift
            X((2 + ((n * (p_i - 1) + 1):(n * p_i))), :) = y(:, (p + 1 - p_i):(T - p_i));
        else
            X((1 + ((n * (p_i - 1) + 1):(n * p_i))), :) = y(:, (p + 1 - p_i):(T - p_i));
        end
    end
    Pi = Y * X' / (X * X');
    e = Y - Pi * X;
    k = (n ^ 2) * p;
end