# EXT-Tseries-Meteoceanic
Simulation of univariate extreme time series accounting for temporal dependence

For better modelling surge-induced coastal flooding, we analyse extreme time series and build a simulator of extreme time series. 
Our main contributions are the following:
- Accounting for temporal dependence and short-tailed behavior. Use of regularly varying functions, extreme value model and polar decomposition. 
- Dimensionality reduction to simulate new extreme time series based on copula models. 
- Tunable aspects, allowing to produce consecutive extremes.
- Several methods proposed to validate the simulation method, using extreme value theory, PCA decomposition and classification two-samples test. 

To run the codes, first run [Detrending_Season](./Detrending_Season.Rmd) for detrending the time series (S) then run [Whitening](./Whitening.Rmd) for whitening them. Finally, using the residuals (./residuals/Name_variable_residuals.csv), run [Simul_ext_Residuals] (./Simul_ext_Residuals.Rmd) and [Simul_Ext_Surge](/Simul_Ext_Surge..Rmd) for simulating respectively extreme residuals and extreme time series of interest. 

All details of the methods are provided in [Gorse et al. (2025)](ajout du lien vers le preprint!!). The methods are applied to [surge data](./Data_Surge.csv) of Gavres site (French Atlantic coast). If you use this dataset, please refer to [Idier et al. (2020)](https://link.springer.com/article/10.1007/s11069-020-03882-4).
