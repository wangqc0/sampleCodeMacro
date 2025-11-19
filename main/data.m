% This is a replicable code sample of downloading US economic data from
% FRED database.
clear

% initialization
cd ..
addpath('function')
if ~isfolder('cleaned')
    mkdir('cleaned')
end

%% Download data

GDPC1 = readtable(['https://fred.stlouisfed.org/graph/fredgraph.csv?' ...
    'id=GDPC1&cosd=1971-01-01&coed=2024-12-31'], 'FileType', 'text');
CPIAUCSL = readtable(['https://fred.stlouisfed.org/graph/fredgraph.csv?' ...
    'id=CPIAUCSL&cosd=1971-01-01&coed=2024-12-31'], 'FileType', 'text');
FEDFUNDS = readtable(['https://fred.stlouisfed.org/graph/fredgraph.csv?' ...
    'id=FEDFUNDS&cosd=1971-01-01&coed=2024-12-31'], 'FileType', 'text');
TWEXBMTH = readtable(['https://fred.stlouisfed.org/graph/fredgraph.csv?' ...
    'id=TWEXBMTH&cosd=1971-01-01&coed=2006-01-01'], 'FileType', 'text');
TWEXBGSMTH = readtable(['https://fred.stlouisfed.org/graph/fredgraph.csv?' ...
    'id=TWEXBGSMTH&cosd=2006-01-01&coed=2024-12-31'], 'FileType', 'text');

%% Adjust data

% GDP
GDPC1_quarterly = table2timetable(GDPC1);
% CPI
CPIAUCSL_quarterly = convert2quarterly(table2timetable(CPIAUCSL), 'Aggregation', 'mean');
CPIAUCSL_quarterly.observation_date = dateshift(CPIAUCSL_quarterly.observation_date, 'start', 'quarter');
% FFR
FEDFUNDS_quarterly = convert2quarterly(table2timetable(FEDFUNDS), 'Aggregation', 'mean');
FEDFUNDS_quarterly.observation_date = dateshift(FEDFUNDS_quarterly.observation_date, 'start', 'quarter');
% EX
convert_ratio_TWEXBGSMTH_to_TWEXBMTH = TWEXBMTH.TWEXBMTH(TWEXBMTH.observation_date == '2006-01-01') / TWEXBGSMTH.TWEXBGSMTH(TWEXBGSMTH.observation_date == '2006-01-01');
EX = synchronize(table2timetable(TWEXBMTH), table2timetable(TWEXBGSMTH));
EX.EX_converted = EX.TWEXBMTH;
EX.EX_converted(EX.observation_date >= "2006-01-01") = EX.TWEXBGSMTH(EX.observation_date >= "2006-01-01") * convert_ratio_TWEXBGSMTH_to_TWEXBMTH;
EX.EX_converted_sa = seasonAdj(EX.EX_converted, 12, [3 5], 'multiplicative');
EX_quarterly = convert2quarterly(EX(:, "EX_converted_sa"), 'Aggregation', 'mean');
EX_quarterly.observation_date = dateshift(EX_quarterly.observation_date, 'start', 'quarter');
% Aggregate
economy_us = synchronize(GDPC1_quarterly, CPIAUCSL_quarterly);
economy_us = synchronize(economy_us, FEDFUNDS_quarterly);
economy_us = synchronize(economy_us, EX_quarterly);

%% Save data

save('cleaned/economy_us.mat', 'economy_us')
