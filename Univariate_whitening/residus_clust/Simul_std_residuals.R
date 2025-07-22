rm(list=ls())
set.seed(133)
par(mfrow=c(1,1))

source("fonctions.R")
source("fonctions_Pareto.R")
require(reshape2)
require(ggplot2)
require(gridExtra)
require(grid)
require(dplyr)
require(threshr)
require(extRemes)
require(ggside)

cols_<-c("data"="blue","simul"="orange","confidence_band"="darkblue")

# Obtenir à partir des simulations  ---------------------------------------
# des résidus de nouvelles trajectoires -----------------------------------
#################

Nom_v<-"Surcote"
inds_extremes<-paste0("residus_clust/inds_extremes_donnees_",Nom_v,".csv")
Individus_exts<-read.csv(inds_extremes)[,2]
lien_DF<-paste0("Travail_hypothèse_RV/EVA_",Nom_v,"_residus.csv")
DF<-read.csv(lien_DF)

# Pareto ------------------------------------------------------------------
Pareto_obs<-read.csv(file=paste0("residus_clust/",Nom_v,"_obs_ech_Pareto.csv"))
df<-Pareto_obs[,c(2:38)]
NL2_P<-apply(df[Individus_exts,],MARGIN = 1,FUN = calcul_norme_L2)
seuil_P<-min(NL2_P)

# Extract some trajectories -----------------------------------------------
M<-1000

elt<-sapply(df, unlist)
Sampling<-sample(c(1:nrow(elt)),size=M)
df_ss_base<-elt[Sampling,]

Locations<-matrix(c(1:37)/37,nrow = 37,ncol = 1,byrow = FALSE)
Madogram<-SpatialExtremes::fmadogram(data=df_ss_base,Locations)

modelisation<-SpatialExtremes::fitmaxstab(data = df_ss_base,Locations,cov.mod="brown")
modelisation$param

# Hs-Smooth proche de 2 (même chose que stage) ----------------------------

# Surcote- ----------------------------------------------------------------

par(mfrow=c(1,1))
values<-sapply(Locations,function(x){return(modelisation$ext.coeff(x))})
plot(Locations,values)

Preds_BR<-fonction_simul_BR_accept_reject(NB_BR = 1000,modele_BR = modelisation,
                                          seuil_Frech = seuil_P)

inds_ss_base<-which(apply(X = df_ss_base,MARGIN = 1,FUN = calcul_norme_L2)>seuil_P)
par(mfrow=c(1,2))
matplot(t(df_ss_base[-inds_ss_base,]),type="l")
matplot(t(Preds_BR),type="l")
par(mfrow=c(1,1))

# Conversion --------------------------------------------------------------
vect_AD<-rep(NA,ncol(Preds_BR))
for(t in c(1:ncol(Preds_BR))){
  vect_AD[t]<-goftest::ad.test(x = Preds_BR[,t],null = extRemes::"pevd",
                   type="GEV",shape=1,loc=1,scale=1)$p.value
}
plot(vect_AD)

# Rejet hypothese nulle -----------------------------------------------------
Sans_condition<-SpatialExtremes::rmaxstab(1000,Locations,cov.mod="brown",
                          range=modelisation$param[1], 
                          smooth=modelisation$param[2])
vect_AD2<-rep(NA,ncol(Preds_BR))
for(t in c(1:ncol(Preds_BR))){
  vect_AD2[t]<-goftest::ad.test(x = Sans_condition[,t],null = extRemes::"pevd",
                               type="GEV",shape=1,loc=1,scale=1)$p.value
}
plot(vect_AD2)
# Conserve hypothese nulle ------------------------------------------------


# conversion_uniforme -----------------------------------------------------
Unifs_BR<-exp(-Preds_BR^(-1))

# conversion de uniforme a orig -------------------------------------------
DF<-read.csv(file=paste0("Travail_hypothèse_RV/EVA_",Nom_v,"_residus.csv"))
Residus_donnees<-read.csv(file=paste0("residus_clust/",Nom_v,"_residus.csv"))[,c(2:38)]
colnames(Residus_donnees)<-c(1:ncol(Residus_donnees))
NDENS<-2^(14)
Orig_BR<-matrix(NA,nrow = nrow(Unifs_BR),
                 ncol = ncol(Unifs_BR))
for(t in c(1:ncol(Residus_donnees))){
  p_u<-DF$p_u_t[t]
  u_f<-DF$seuil_t[t]
  variable<-Residus_donnees[,t]
  kernel_DENS<-density(x=variable,n=NDENS)
  step<-diff(kernel_DENS$x)[1]
  integrale<-cumsum(kernel_DENS$y)*step
  val_orig<-Unifs_BR[,t]
  inds_exts<-which(Unifs_BR[,t]>(1-p_u))
  val_orig[inds_exts]<-qgp_craft((Unifs_BR[inds_exts,t]-1+p_u)/p_u,sigma = DF$forme_t[t],
                                xi = DF$forme_t[t])+u_f
  Indices_d_min<-sapply(Unifs_BR[-inds_exts,t],function(x){
    return(which.max(x-integrale<0)-1)})
  val_orig[-inds_exts]<-as.numeric(kernel_DENS$x[Indices_d_min])
  Orig_BR[,t]<-val_orig
}

par(mfrow=c(1,2))
matplot(t(Orig_BR),type="l")
matplot(t(Residus_donnees[-Individus_exts,]),type="l")
par(mfrow=c(1,1))

# Projections -------------------------------------------------------------
Simul_Epsi_std<-Orig_BR
Real_Epsi_std<-Residus_donnees[-Individus_exts,]
mu_std<-colMeans(x = Real_Epsi_std)
sd_std<-apply(X = Real_Epsi_std,MARGIN = 2,FUN = sd)
PCA_std<-FactoMineR::PCA(X =Real_Epsi_std,scale.unit = TRUE,
                         graph = TRUE)
plot(1-c(0,cumsum(PCA_std$eig[,1])/sum(PCA_std$eig[,1])))
Simul_Epsi_std_for_PCA<-scale(Simul_Epsi_std,center = mu_std,
                              scale = sd_std)
Coords_simul_sdt<-Simul_Epsi_std_for_PCA%*%PCA_std$svd$V[,c(1:2)]
Coords_real<-PCA_std$ind$coord[,c(1:2)]
ks.test(Coords_real[,1],Coords_simul_sdt[,1])

# Conserver hypothese nulle pour Surcote ------------------------------------
# Pas le cas pr Surcote. --------------------------------------------------

