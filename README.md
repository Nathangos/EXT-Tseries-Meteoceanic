# EXT-Tseries-Meteoceanic
Simulation of univariate extreme time series accounting for temporal dependence

For better modelling surge-induced coastal flooding, we analyse extreme time series and build a simulator of extreme time series. 
Our main contributions are the following:
- Accounting for temporal dependence and short-tailed behavior. Use of regularly varying functions, extreme value model and polar decomposition. 
- Dimensionality reduction to simulate new extreme time series based on copula models. 
- Tunable aspects, allowing to produce consecutive extremes.
- Several methods proposed to validate the simulation method, using extreme value theory, PCA decomposition and classification two-samples test. 

To run the codes, first run [Detrending_Season](./Detrending_Season.Rmd) for detrending the time series then run [Whitening](./Whitening.Rmd) for whitening them.
Finally, using the residuals displayed in folder (./residuals), run [Simul_ext_Residuals](./Simul_ext_Residuals.Rmd) and [Simul_Ext_Forcing](./Simul_Ext_Forcing.Rmd) for simulating respectively extreme residuals and extreme time series of interest. Additional results are displayed in [Complements_analysis](./Complements_analysis.Rmd). 

All details of the methods are provided in [Gorse et al. (2025)](https://arxiv.org/abs/2508.13687). The methods are applied to [surge data](./Data/Data_Surge.csv) and [wind_speed data](./Data/Data_U.csv) of Gavres site (French Atlantic coast). If you use these datasets, please refer to [Idier et al. (2020)](https://link.springer.com/article/10.1007/s11069-020-03882-4).
