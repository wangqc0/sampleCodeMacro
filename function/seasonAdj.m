function dt = seasonAdj(y, s, m, decomp)
%SEASONADJ Seasonal adjustment for time series data, which resembles the
% process used within the X-12-ARIMA seasonal adjustment program of the US
% Census Bureau
% Reference: https://www.mathworks.com/help/econ/seasonal-adjustment-1.html
% Reference: Ladiray and Quenneville (2001), Seasonal Adjustment with the
% X-11 Method, ISBN: 9780387951713
% Reference: Doherty (2001), The Surrogate Henderson Filters in X-11,
% Australian and New Zealand Journal of Statistics 43(4), 2001, 385–392
%   y: input timeseries
%   s: periodicity (12 for monthly data with annual periodicity; s can be
%    4, 6, 8, 12, 22)
%   m: parameter in the seasonal filter S(n, m) for the detrended series
%    (n = 3; m can be 1, 3, 5, 9; stable filter when m = 1)
%   decomp: decomposition method ('additive' or 'multiplicative'; additive
%    by default)
    s_available = [4; 6; 8; 12; 22];
    if ~any(s_available == s)
        error(['Inconsistent input for periodicity: should be 4, 6, 8, 12' ...
            ' or 22'])
    end
    if length(m) ~= 2
        error(['Inconsistent input for seasonal filter parameter:' ...
            'should be 1 for stable filter and 3, 5 or 9 for MA filter'])
    end
    if ~exist('decomp', 'var')
        decompfun = @(x, y)(x - y);
    elseif strcmp(decomp, 'additive')
        decompfun = @(x, y)(x - y);
    elseif strcmp(decomp, 'multiplicative')
        decompfun = @(x, y)(x ./ y);
    else
        error(['Inconsistent input for decomposition method:' ...
            'should be either additive or multiplicative'])
    end
    T = length(y);
    N = s + 1;
    % detrend the data using a N-term moving average
    % 1. Obtain a first estimate of the trend component
    sWN = [1 / (s * 2); repmat(1 / s, s - 1, 1) ; 1 / (s * 2)];
    yS = conv(y, sWN, 'same');
    yS(1:(s / 2)) = yS(s/2 + 1);
    yS((T - s / 2 + 1):T) = yS(T - s / 2);
    % 2. Detrend the original series
    xt = decompfun(y, yS);
    % 3. Apply a seasonal filter to the detrended series
    % create seasonal indices
    sidx = cell(s,1);
    for i = 1:s
        sidx{i, 1} = i:s:T;
    end
    % apply filter
    if m(1) == 1
        sst = cellfun(@(x) mean(xt(x)), sidx);
        % put smoothed values back into a vector of length N
        nc = floor(T / s);
        rm = mod(T, s);
        sst = [repmat(sst, nc, 1); sst(1:rm)];
        % center the seasonal estimate
        sBar = mean(sst);
        s1 = decompfun(sst, sBar);
    else
        switch m(1)
            case 3
                sW1 = [1; 2; 3; 2; 1] / 9;
                aW1 = [7 11; 10 11; 7 5; 3 0] / 27;
            case 5
                sW1 = [1; 2; 3; 3; 3; 2; 1] / 15;
                aW1 = [9 15 17; 13 15 17; 13 15 17; 13 11 9; 8 4 0; 4 0 0] / 60;
            case 9
                sW1 = [1; 2; repmat(3, 7, 1); 2; 1] / 27;
                aW1 = [86 145 177 213 252; 123 141 167 197 227; ...
                    121 135 158 181 202; 120 131 147 164 177; ...
                    119 126 136 148 115; 117 120 136 94 52; ...
                    116 116 81 29 0; 114 77 33 0 0; ...
                    75 35 0 0 0; 35 0 0 0 0] / 1026;
        end
        % apply filter to each month
        shat = NaN * y;
        for i = 1:s
            ns = length(sidx{i});
            first = 1:(m(1) + 1);
            last = (ns - m(1)):ns;
            dat = xt(sidx{i});
            sd = conv(dat, sW1, 'same');
            sd(1:((m(1) + 1) / 2)) = conv2(dat(first), 1, rot90(aW1, 2), 'valid');
            sd((ns - (m(1) - 1) / 2):ns) = conv2(dat(last), 1, aW1, 'valid');
            shat(sidx{i}) = sd;
        end
        % N-term moving average of filtered series
        sb = conv(shat, sWN, 'same');
        sb(1:(s / 2)) = sb((s + 1):(s + (s / 2)));
        sb((T - s / 2 + 1):T) = sb((T - s / 2 + 1 - s):(T - s));
        % 4. Deseasonalize the original series
        % center to get final estimate
        s1 = decompfun(shat, sb);
    end
    % 5. Obtain a second estimate of the trend component
    % apply N-term Henderson filter
    dt = decompfun(y, s1);
    p = s / 2;
    n = p + 2;
    i = [flip(1:p) (0:p)]';
    sWH = 315 * ((n - 1)^2 - i.^2) .* (n^2 - i.^2) .*  ((n + 1)^2 - i.^2) .* (3*n^2 - 16 - 11 * i.^2) ...
        / (8 * n * (n^2 - 1) * (4 * n^2 - 1) * (4 * n^2 - 9) * (4 * n^2 - 25));
    aWH = zeros(s, s / 2);
    switch N
        case 5
            R = .001;
        case 7
            R = 4.5;
        case 9
            R = 1;
        case 13
            R = 3.5;
        case 23
            R = 4.5;
    end
    D = 4 / pi / R ^ 2;
    for j = 1:(s / 2)
        M = N - j;
        aWH(1:M, j) = flip(sWH(1:M) + sum(sWH((M + 1):N)) / M + ...
            ((1:M)' - (M + 1) / 2) * D / (1 + M * (M - 1) * (M + 1) / 12 * D) * sum((((M + 1):N) - (M + 1) / 2) * sWH((M + 1):N)));
    end
    % apply N-term Henderson filter
    first = 1:s;
    last = (T - s + 1):T;
    hN = conv(dt, sWH, 'same');
    hN((T - (s / 2) + 1):end) = conv2(dt(last), 1, aWH, 'valid');
    hN(1:(s / 2)) = conv2(dt(first), 1, rot90(aWH, 2), 'valid');
    % hN is the final estimate of the trend component
    % 6. Detrend the original series again
    % new detrended series
    xt = decompfun(y, hN);
    % 7. Apply a seasonal filter to the detrended series
    % apply filter again
    if m(2) == 1
        sst = cellfun(@(x) mean(xt(x)), sidx);
        % put smoothed values back into a vector of length N
        nc = floor(T / s);
        rm = mod(T, s);
        sst = [repmat(sst, nc, 1); sst(1:rm)];
        % center the seasonal estimate
        sBar = mean(sst);
        s2 = decompfun(sst, sBar);
    else
        switch m(2)
            case 3
                sW2 = [1; 2; 3; 2; 1] / 9;
                aW2 = [7 11; 10 11; 7 5; 3 0] / 27;
            case 5
                sW2 = [1; 2; 3; 3; 3; 2; 1] / 15;
                aW2 = [9 15 17; 13 15 17; 13 15 17; 13 11 9; 8 4 0; 4 0 0] / 60;
            case 9
                sW2 = [1; 2; repmat(3, 7, 1); 2; 1] / 27;
                aW2 = [86 145 177 213 252; 123 141 167 197 227; ...
                    121 135 158 181 202; 120 131 147 164 177; ...
                    119 126 136 148 115; 117 120 136 94 52; ...
                    116 116 81 29 0; 114 77 33 0 0; ...
                    75 35 0 0 0; 35 0 0 0 0] / 1026;
        end
        % apply filter to each month
        shat = NaN * y;
        for i = 1:s
            ns = length(sidx{i});
            first = 1:(m(2) + 1);
            last = (ns - m(2)):ns;
            dat = xt(sidx{i});
            sd = conv(dat, sW2, 'same');
            sd(1:((m(2) + 1) / 2)) = conv2(dat(first), 1, rot90(aW2, 2), 'valid');
            sd((ns - (m(2) - 1) / 2):ns) = conv2(dat(last), 1, aW2, 'valid');
            shat(sidx{i}) = sd;
        end
        % 8. Deseasonalize the original series
        % N-term moving average of filtered series
        sb = conv(shat, sWN, 'same');
        sb(1:(s / 2)) = sb((s + 1):(s + (s / 2)));
        sb((T - s / 2 + 1):T) = sb((T - s / 2 + 1 - s):(T - s));
        % center to get final estimate
        s2 = decompfun(shat, sb);
    end
    % s2 is the final estimate of the seasonal component
    % decentralized series
    dt = decompfun(y, s2);
end
