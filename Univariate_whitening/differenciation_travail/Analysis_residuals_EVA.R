rm(list=ls())
set.seed(133)

#### Packages. #####
require(corrplot)
require(ggplot2)
require(dplyr)
require(reshape2)
require(FactoMineR)
require(fda)
require(extRemes)
require(parallel)
require(POT)
require(tibble)
require(threshr)
require(tseries)

### Source des fonctions. ####
source("fonctions/fonctions.R")
F_names<-ls()
nb_coeurs<-detectCores()-4
coeurs<-makeCluster(nb_coeurs)
clusterExport(coeurs,varlist=F_names)


nom_variable<-"Surcote"
Residus<-read.csv(file = paste0("residus_clust/",nom_variable,
                                "_residus.csv"))[,c(2:38)]

# EVA L2 norm -------------------------------------------------------------
L2_norm_resid<-apply(X = Residus,MARGIN = 1,FUN = calcul_norme_L2)
vect_k<-c(20:600)
Nom_graph<-ifelse(nom_variable=="Surcote","S",
                  nom_variable)
# TGraphique_resid<-paste0("the L2 norm of ", 
#                          Nom_graph," (residual)")
# Graphiques_estimateurs_gamma(series = L2_norm_resid,
#                              vect_k = vect_k,
#                              Titre_graphique = NULL,
#                              NB_annees = 37)

vect_temps<-c(19)
C_TH<--0.5

Fct_test_ext_Reiss<-function(lag,series,C_th){
  fin<-length(series)-lag
  x_<-series[c(1:fin)]
  debut<-lag+1
  x_lag1<-series[c(debut:length(series))]
  #Reiss
  Test_lag1<-extRemes::taildep.test(x_,x_lag1,cthresh = C_th)
  u<-do.call("relative.rank",c(list(x=x_)))-1
  v<-do.call("relative.rank",c(list(x=x_lag1)))-1
  ci<-u+v
  ci<-ci[ci>C_th]
  resultat_base<-ci/C_th

  # beta estimator ----------------------------------------------------------
  estim_beta_Moment<-(1-mean(resultat_base))^(-1)-2
  
  # Anderson Darling test ---------------------------------------------------
  Unif<-resultat_base^(estim_beta_Moment+1)
  AD_pval<-goftest::ad.test(x = Unif,null = "punif")$p.value

  return(list("pval"=Test_lag1$p.value,
         "sample"=resultat_base,
         "ci"=ci,"AD_beta"=AD_pval,
         "beta"=estim_beta_Moment))
}
Mat_ind_ext_p_valeur_Reiss<-list()
for(col in vect_temps){
  
  resids_col<-Residus[,col]
  L<-length(resids_col)-1
  df<-cbind.data.frame(resids_col[1:L],
                       resids_col[2:length(resids_col)])
  par(mfrow=c(1,2))
  # chi measure -------------------------------------------------------------
  POT::chimeas(df,which = 1)
  abline(h=0,col="red")
  POT::chimeas(df,which = 2)
  abline(h=0,col="red")
  nom_graph<-ifelse(test = nom_variable=="Surcote",
                    "S",nom_variable)
  # Titre_chi<-paste0("Analyis at time t=",col,
  #                   " for ",nom_graph, " (lag=1)")
  # mtext(text =Titre_chi,outer = TRUE,line=-3)
  par(mfrow=c(1,1))

  # Indep tail --------------------------------------------------------------
  # Depends on an exceedance rate -------------------------------------------

  Mat_ind_ext_p_valeur_Reiss[[paste0("lag_",col)]]<-sapply(c(1:50),Fct_test_ext_Reiss,
                                                           series=resids_col,C_th=C_TH)
}
J<-ncol(Mat_ind_ext_p_valeur_Reiss$lag_1)
vect_beta<-rep(NA,J)
vect_Ad<-rep(NA,J)
for(j in c(1:J)){
  vect_beta[j]<-Mat_ind_ext_p_valeur_Reiss$lag_1[,j]$beta
  vect_Ad[j]<-Mat_ind_ext_p_valeur_Reiss$lag_1[,j]$AD_beta
}
plot(vect_Ad,vect_beta,main="beta estimator for various lags",
     ylab=expression(beta),xlab="p val AD test")

J<-ncol(Mat_ind_ext_p_valeur_Reiss$lag_1)
vect_beta<-rep(NA,J)
vect_Ad<-rep(NA,J)
for(j in c(1:J)){
  vect_beta[j]<-Mat_ind_ext_p_valeur_Reiss$lag_7[,j]$beta
  vect_Ad[j]<-Mat_ind_ext_p_valeur_Reiss$lag_7[,j]$AD_beta
}
plot(vect_Ad,vect_beta,main="beta estimator for various lags",
     ylab=expression(beta),xlab="p val AD test")

# Fin code ----------------------------------------------------------------
# ------------------------------------------------------------------------

