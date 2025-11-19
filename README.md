# Sample Code in Macroeconomic Analysis

By: Qichao Wang

Last update on: 19 November 2025 (Version 1)

A sample project containing downloading macroeconomic data from FRED and applying VAR and FEVD analysis using MATLAB.

## Project Structure

```
Qichao-Wang-code-sample/
├── main/                       # MATLAB scripts
│   ├── run_all.m               # Run all scripts
│   ├── data.m                  # Download data from FRED
│   └── analysis.m              # VAR and FEVD analysis
└── function/                   # MATLAB functions
    ├── seasonAdj.m             # Seasonal adjustment on time series
    ├── var_ols.m               # VAR: time series -> estimate
    ├── var_irf_chol.m          # VAR: estimate -> IRF value
    ├── var_irf_chol_bs.m       # VAR: estimate -> booststrapped IRF draws
    ├── var_irf_plot.m          # VAR: IRF value + booststrapped IRF draws -> IRF plot
    ├── fevd.m                  # FEVD: VAR estimate -> FEVD value
    ├── fevd_plot.m             # FEVD: FEVD value -> FEVD plot
    ├── num_sprintf.m           # Numeric value -> Formatted string
    ├── latexTable.m            # Table class object -> LaTex string cell (modified based on Duenisch, 2016)
    ├── latexTable_add_row_bottom.m
    └── latexTable_save.m
```

## Features

### Implemented Methods

1. **Seasonal Adjustment** (MATLAB)
   - Simple implementation of the X-12-ARIMA seasonal adjustment program of the US Census Bureau

2. **Vector Autoregression (VAR)** (MATLAB)
   - Cholesky decomposition to obtain the variance matrix
   - Use matrix operation to speed up relative to for loop when bootstrapping
   - Adjust the error band based on the difference between the point estimate and the bootstrap median

3. **Forecast Error Variance Decomposition (FEVD)** (MATLAB)
   - Use matrix operations to speed up relative to a for loop
   - Customizable color on the plot
   - Adjust the error band based on the difference between the point estimate and the drawn median

### Datasets

The data used in the project is downloaded from the [FRED](https://fred.stlouisfed.org) database provided by the Federal Reserve Bank of St. Louis.

## Usage

### Running Analyses

#### MATLAB Analyses

```matlab
# Run all
run_all
```

The script `run_all.m` automatically and sequentially runs the data downloading and cleaning script (`data.m`) and then the analysis script (`analysis.m`). It generates the following files saved in the respectively generated folders:
1. Cleaned data in `.mat` format (saves to `cleaned/`)
2. Figures in `.eps` format (saves to `figure/`)
3. Tables in `.tex` format (saves to `table/`)

## License

This is a research project for educational and training purposes.
