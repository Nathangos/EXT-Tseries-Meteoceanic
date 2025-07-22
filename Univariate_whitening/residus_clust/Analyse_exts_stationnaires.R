rm(list=ls())
set.seed(133)
par(mfrow=c(1,1))
source("fonctions.R")
require(reshape2)
require(ggplot2)
require(gridExtra)
require(grid)
require(dplyr)
require(extRemes)

nom_variable<-"Hs"
repertoire<-"ss_tend/"
type_donnees<-"HIVER"
if(type_donnees=="HIVER"){
  repertoire<-paste0(repertoire,"HIVER/")
}
if(type_donnees=="HORS_HIVER"){
  repertoire<-paste0(repertoire,"HORS_HIVER/")
}
nom_recherche<-paste0(repertoire,nom_variable,"_ss_tend.csv")
obs<-read.csv(nom_recherche)[,2:38]
colnames(obs)<-c(1:37)
borne_passe<-nrow(obs)-1
obs_passe<-obs[c(1:borne_passe),]
obs_present<-obs[c(2:nrow(obs)),]

# modele_ar ---------------------------------------------------------------
fich_reg<-paste0("differenciation_travail/",nom_variable,"_Reg_AR1.csv")
Mod_AR<-read.csv(file=fich_reg)
Fin<-ncol(Mod_AR)
Mod_AR<-Mod_AR[,c(2:Fin)]
j<-19
Fin<-ncol(Mod_AR)-1
vect_regression<-Mod_AR[j,c(1:Fin)]
Constante<-c(unlist(Mod_AR$Intercept)[j]%*%(1-sum(unlist(vect_regression))))

if(length(vect_regression)==1){
  objet<-diag(obs_passe[,j])
  vecteur_beta1<-rep(vect_regression,length(obs_passe[,j]))
  resultat_AR<-c(t(objet%*%vecteur_beta1))
  X_projection<-resultat_AR+rep(Constante,length(obs_passe))
}
plot(X_projection,obs_present[,j])

Residus<-read.csv(paste0("residus_clust/",nom_variable,"_residus.csv"))[,2:38]
colnames(Residus)<-c(1:ncol(Residus))
# Extremes analysis -------------------------------------------------------
P_u<-0.10
series_j<-obs_present[,j]
seuil_t<-as.numeric(quantile(series_j,1-P_u))

# sans_cov ----------------------------------------------------------------
series_j_ext<-subset(series_j,series_j>seuil_t)

# test stationnarite ------------------------------------------------------
tseries::adf.test(series_j_ext)
tseries::kpss.test(series_j_ext)
acf(series_j_ext)
pacf(series_j_ext)

resultat_log<-extRemes::fevd(x = series_j,threshold = seuil_t, 
               type = "GP")
plot(resultat_log)
objet<-summary(resultat_log)
BIC<-objet$BIC
BIC

# avec_cov-projection ----------------------------------------------------------------
data_EVD<-cbind.data.frame(series_j,(X_projection-mean(X_projection))/sd(X_projection))
colnames(data_EVD)<-c("value_present","predictions_k_passe")
NPY<-round(length(series_j)/37)
resultat_log_covariable<-extRemes::fevd(value_present,data_EVD,threshold = seuil_t, 
                                        threshold.fun = ~predictions_k_passe,
                                        type = "GP",units = paste0(NPY,"/year"))
plot(resultat_log_covariable)

resultat_log_covariable$results$par
objet_cov<-summary(resultat_log_covariable)
BIC_cov<-objet_cov$BIC
BIC_cov
ci(resultat_log_covariable,type="parameter")


# avec cov-residual -------------------------------------------------------
resid_j<-Residus[c(2:nrow(obs)),j]
data_EVD_residual<-cbind.data.frame(series_j,(resid_j-mean(resid_j))/sd(resid_j))
colnames(data_EVD_residual)<-c("value_present","residuals_AR")

NPY<-round(length(series_j)/37)
resultat_log_covariable<-extRemes::fevd(value_present,data_EVD_residual,threshold = seuil_t, 
                                        threshold.fun = ~residuals_AR,
                                        scale.fun = ~residuals_AR,
                                        type = "GP",units = paste0(NPY,"/year"))
plot(resultat_log_covariable)
BIC_cov_resid<-summary(resultat_log_covariable)$BIC
BIC_cov_resid
