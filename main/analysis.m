% This is a replicable code sample of VAR and FEVD methods using input
% empirical data.
% The output include figures and tables in the corresponding subfolders.
clear

% initialization
cd ..
addpath('cleaned')
addpath('function')
if ~isfolder("figure")
    mkdir("figure")
end
if ~isfolder("table")
    mkdir("table")
end

% parameter
M = 4;
Q = 4;
p_min = 4;
p_max = 15;
h_max = 20;
n_mc = 1000;
n_bs = 1000;
T_burnin = 50;
irf_band_confidence = 90;

% set seed
rng(0, 'twister')

% read data
load('cleaned/economy_us.mat')

%% Generate variables

economy_us.gdp = log(economy_us.GDPC1);
economy_us.cpi = log(economy_us.CPIAUCSL);
economy_us.ex = log(economy_us.EX_converted_sa);
economy_us.ffr = economy_us.FEDFUNDS;
economy_us.dgdp = [nan; 100 * 4 * diff(economy_us.gdp)];
economy_us.dcpi = [nan; 100 * 4 * diff(economy_us.cpi)];
economy_us.dex = [nan; 100 * 4 * diff(economy_us.ex)];

%% VAR

% Time series plot
var_select = {'dgdp', 'dcpi', 'dex', 'ffr'};
gdpcpiexffr = economy_us(~any(isnan(economy_us{:, var_select}), 2), var_select);
stackedplot(gdpcpiexffr, 'Color', 'black')
title('GDP, CPI, EX, and FFR (Annualized %Change)')
saveas(gcf, 'figure/figure_1.eps')
close;

% Run VAR to find the optimal number of lags
T = size(gdpcpiexffr, 1) - p_max;
n = size(gdpcpiexffr, 2);
y = gdpcpiexffr{:, :}';
aic_p = nan(p_max - p_min + 1, 1);
for p = p_min:p_max
    [~, e, k] = var_ols(y(:, (1 + p_max - p):end), p);
    var_e = cov(e');
    aic = log(det(var_e)) + 2 * k / T;
    aic_p(p - p_min + 1) = aic;
end
p = find(aic_p == min(aic_p)) + p_min - 1;
% Obtain the residual from VAR
[n, T] = size(y);
[Pi, e, k] = var_ols(y, p);
y_hat = y(:, (p + 1):end) - e;
var_e = cov(e');

% Obtain error band of VAR through Cholesky decomposition
irf_varname = {'GDP (cumulative)', 'CPI (cumulative)', 'EX (cumulative)', 'FFR'};
[dy_dv, dy_de, P] = var_irf_chol(Pi, e, h_max);
[dy_dv_bs, dy_de_bs, P_bs] = var_irf_chol_bs(y, Pi, e, h_max, T_burnin, n_bs);
% Productivity shock
ind_shock = 1;
size_shock = 1;
var_irf_plot(dy_dv, ind_shock, size_shock, 1:3, dy_dv_bs, 1, irf_varname, irf_band_confidence)
sgtitle({'IRF of VAR, Productivity shock', ['Bootstrapping (B = ', num2str(n_bs), ') with bias correction, p = ', num2str(p)]})
saveas(gcf, 'figure/figure_2_1.eps')
close;
% Price shock
ind_shock = 2;
size_shock = 1;
var_irf_plot(dy_dv, ind_shock, size_shock, 1:3, dy_dv_bs, 1, irf_varname, irf_band_confidence)
sgtitle({'IRF of VAR, Price shock', ['Bootstrapping (B = ', num2str(n_bs), ') with bias correction, p = ', num2str(p)]})
saveas(gcf, 'figure/figure_2_2.eps')
close;
% Exchange rate shock
ind_shock = 3;
size_shock = 1;
var_irf_plot(dy_dv, ind_shock, size_shock, 1:3, dy_dv_bs, 1, irf_varname, irf_band_confidence)
sgtitle({'IRF of VAR, Exchange rate shock', ['Bootstrapping (B = ', num2str(n_bs), ') with bias correction, p = ', num2str(p)]})
saveas(gcf, 'figure/figure_2_3.eps')
close;
% FFR shock
ind_shock = 4;
size_shock = 1;
var_irf_plot(dy_dv, ind_shock, size_shock, 1:3, dy_dv_bs, 1, irf_varname, irf_band_confidence)
sgtitle({'IRF of VAR, FFR shock', ['Bootstrapping (B = ', num2str(n_bs), ') with bias correction, p = ', num2str(p)]})
saveas(gcf, 'figure/figure_2_4.eps')
close;

% FEVD
fevd_var = fevd(dy_de, P, 1:3);
fevd_plot(fevd_var, irf_varname)
saveas(gcf, 'figure/figure_3.eps', 'epsc')
close;
% Confidence interval of FEVD (Price shock)
fevd_var_bs = nan(h_max + 1, n, n, n_bs);
for i_bs = 1:n_bs
    fevd_var_bs(:, :, :, i_bs) = fevd(dy_de_bs(:, :, :, i_bs), P_bs(:, :, i_bs), 1:3);
end
fevd_var_bs_upper = prctile(fevd_var_bs, 100 - (100 - irf_band_confidence) / 2, 4);
fevd_var_bs_lower = prctile(fevd_var_bs, (100 - irf_band_confidence) / 2, 4);
fevd_var_bs_mean = mean(fevd_var_bs, 4);
ind_shock = 2;
fevd_var_ffr_band = cat(3, fevd_var(:, :, ind_shock), reshape([fevd_var_bs_upper(:, :, ind_shock) fevd_var_bs_lower(:, :, ind_shock)], h_max + 1, n, 2) + fevd_var(:, :, ind_shock) - fevd_var_bs_mean(:, :, ind_shock));
table_fevd_var_ffr_band = pagetranspose(permute(fevd_var_ffr_band((0:5:h_max) + 1, :, :), [1 3 2]));
table_fevd_var_ffr_band_concatenated = {};
table_fevd_var_ffr_band_concatenated_row = {};
for i_var = 1:numel(irf_varname)
    table_fevd_var_ffr_band_concatenated = [table_fevd_var_ffr_band_concatenated; repmat({''}, 1, h_max / 5 + 1); num_sprintf(splitvars(table(table_fevd_var_ffr_band(:, :, i_var))), {'%0.2f'})];
    table_fevd_var_ffr_band_concatenated_row = [table_fevd_var_ffr_band_concatenated_row; ['\textit{', irf_varname{i_var}, '}']; 'Mean'; [num2str(100 - (100 - irf_band_confidence) / 2) '$\%$']; [num2str((100 - irf_band_confidence) / 2) '$\%$']];
end
table_fevd_var_ffr_band_concatenated_column = cellstr(strcat('$h=', string(0:5:h_max), '$'));
clear input
input.data = table_fevd_var_ffr_band_concatenated{:, :};
input.tableColLabels = table_fevd_var_ffr_band_concatenated_column;
input.tableRowLabels = table_fevd_var_ffr_band_concatenated_row;
input.tablePlacement = 'H';
input.tableBorders = 1;
input.booktabs = 1;
input.tableCaption = 'FEVD statistics: CPI, with GDP, EX, FFR';
input.tableLabel = strcat('tab_3');
input_tab = latexTable(input);
input_tab = latexTable_add_row_bottom(input_tab, {'The bootstrapped results are adjusted according to the differences between the point estimates of the VAR and the bootstrapped means.'}, 'threeparttable', true);
latexTable_save(input_tab, 'table/table_3.tex', 'w')
