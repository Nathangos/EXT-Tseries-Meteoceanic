rm(list=ls())
set.seed(133)
par(mfrow=c(1,1))

#### Import des packages. #####
require(dplyr)
require(ggplot2)
require(reshape2)
require(FactoMineR)
require(extRemes)
require(fda)
require(parallel)
require(VineCopula)
require(grid)
require(gridExtra)
require(mvPot)

### Source des fonctions. ####
source("fonctions_Pareto.R")
source("Modele_ss_tend/Travail_angles.R")
source("fonctions.r")
source("fonctions_r_pareto_EVOLTAU.r")
F_names<-ls()
nb_coeurs<-detectCores()-4


# Ouverture des clusters --------------------------------------------------
coeurs<-makeCluster(nb_coeurs)
clusterExport(coeurs,varlist=F_names)

Simulation_ACP_TS_EXT<-function(type_donnees,coeurs,lien_donnees,nom_variable,p_U,p_L2,n.dens,opt_Frech,M,type_entree,type_vario){
  
}
type_entree<-"residus"
nom_variable<-"Hs"
STD_donnees<-TRUE
P_U<-rep(0.10,37)
PL2<-0.05
opt_Frech<-TRUE
type_V<-"pow"

if(type_entree=="residus"){
  l_donnees<-paste0("residus_clust/",nom_variable,"_residus.csv")
  racine_ML<-c("resultat_ML/residus/")
}
if(type_entree=="donnees_brutes"){
  repertoire<-"ss_tend/HIVER/"
  l_donnees<-paste0(repertoire,nom_variable,"_ss_tend.csv")
  racine_ML<-c("resultat_ML/donnees_brutes/")
}
if(type_entree=="clust"){
  l_donnees<-paste0("residus_clust/clust_",nom_variable,".csv")
  
}

# Import EVA --------------------------------------------------------------

lien_DF<-paste0("Travail_hypothèse_RV/EVA_",nom_variable,"_residus.csv")
DF<-read.csv(lien_DF)

# Import donnees ----------------------------------------------------------

Residus_donnees<-read.csv(file=paste0("residus_clust/",nom_variable,"_residus.csv"))[,c(2:38)]
colnames(Residus_donnees)<-c(1:ncol(Residus_donnees))


# Variance pic maree ------------------------------------------------------

t<-19
periods_years<-c(2,5,10,20,50,80,
                 100,120,200,
                 250,300,500,
                 800)
rate_exces<-DF$p_u_t[t]
NB_years_total<-37
series<-Residus_donnees[,t]
npy<-length(series)/NB_years_total
p_proba<-1-(npy*periods_years*rate_exces)^(-1)

seuil<-DF$seuil_t[t]

shape_opitz<-DF$forme_t[t]
scale_opitz<-DF$échelle_t[t]
series_excedents<-subset(series,series>seuil)
resultat<-Fnct_boostrap_param(B = 500,nb_sample = length(series_excedents),
                    GPD_liste = list("shape"=shape_opitz,"scale"=scale_opitz,
                                     "threshr"=seuil), 
                    p=p_proba,alpha=0.01)
intervalle_conf_pic<-as.numeric(quantile(resultat$shape,
                                         probs = c(0.05,0.95)))
plot(density(resultat$shape),
     main="Densité gamma au pic de la marée ")

# Ensemble params autres temps --------------------------------------------
estimateurs<-DF$forme_t
Min<-min(min(estimateurs),
         min(intervalle_conf_pic))
Max<-max(max(estimateurs),
         max(intervalle_conf_pic))
plot(c(1:nrow(DF)),estimateurs,ylim=c(Min,Max))
abline(h=intervalle_conf_pic,col="red")

# L'hypothèse nulle d'un gamma constant se confirme -----------------------
NL2r<-apply(X = Residus_donnees,MARGIN = 1,FUN = calcul_norme_L2)
mu<-quantile(NL2r,0.95)
indices_ext_scale_orig<-which(NL2r>mu)

# Transformation en Frechet -----------------------------------------------
gpd_p_val<-rep(NA,37)
frech_p_val<-gpd_p_val
indices_finaux<-c(1:nrow(Residus_donnees))
Residus_modif<-Residus_donnees
for(j in c(1:ncol(Residus_modif))){
  series_j<-Residus_donnees[,j]
  seuil<-DF$seuil_t[j]
  echelle<-DF$échelle_t[j]
  indices<-which(series_j>seuil)
  series_GPD<-series_j[indices]
  test_gpd<-goftest::ad.test(x = series_GPD, 
                   null = extRemes::"pevd", 
                   threshold=seuil, 
                   scale=echelle, 
                   shape=shape_opitz, 
                   type="GP")
  gpd_p_val[j]<-test_gpd$p.value
  coeff<-shape_opitz/echelle
  inverse_sh<-shape_opitz^(-1)
  series_F<-(-log(1-(1+(coeff*(series_j-seuil)))^(-inverse_sh)))^(-1)
  Residus_modif[,j]<-series_F
  test<-goftest::ad.test(x = series_F[indices], 
                   null = extRemes::"pevd", 
                   shape=1, 
                   scale=1, 
                   loc=1, 
                   type="GEV")     
  frech_p_val[j]<-test$p.value
  indices_finaux<-intersect(indices_finaux,indices)
}
plot(c(1:37),gpd_p_val)
plot(c(1:37),frech_p_val)  
indices_finaux<-intersect(indices_finaux,indices_ext_scale_orig)
Obs_TX<-Residus_modif[indices_finaux,]
indices_nonext<-which(!c(1:nrow(Residus_modif))%in%indices_ext_scale_orig==TRUE)
print(length(indices_nonext))
L2_tf<-apply(Residus_modif, MARGIN = 1,
                        FUN = calcul_norme_L2)
max(NL2r[indices_nonext])
print(min(NL2r[indices_finaux]))


# bien au-dessus ----------------------------------------------------------

# MAIS  -------------------------------------------------------------------
v_seuil<-max(L2_tf[indices_nonext])
print(v_seuil)
matplot(t(Residus_modif[indices_ext_scale_orig,]),type="l")
matplot(t(Residus_modif[indices_nonext,]),type="l")

# travail avec vario ------------------------------------------------------
LOC<-expand.grid(1:37/37)
vecteur_temps<-c(1:37)/37
Excedents<-as.list(as.data.frame(t(Obs_TX)))

# BResnick ----------------------------------------------------------------
weigthFun <- function(x, u){
  return(x * (1 - exp(-((calcul_norme_L2(x) / u) - 1))))
  #return(x * (1 - exp(-(sum(x) / u) - 1)))
  
}
#Define partial derivative of weighting function
dWeigthFun <- function(x, u){
  return((1 - exp(-((calcul_norme_L2(x) / u) - 1))) + 
           (1/length(x))*((calcul_norme_L2(x))^(-1))*(x^(2)/u)*exp( - ((calcul_norme_L2(x) / u) - 1)))
  #(1 - exp(-(sum(x / u) - 1))) + (x / u) * exp( - (sum(x / u) - 1))
  
}


f_h_param<-function(h,parametre){
  return(parametre[2]*(abs(h))^(parametre[1]))
}
objectiveFunction = function(parameter, exceedances, loc, vario, weigthFun, dWeigthFun, threshold){
  if(parameter[1]<0|parameter[2]<0|parameter[1]>2){return(1e50)}
  varioModel <- function(h){
    vario(h = h,param = parameter)
  }
  mvPot::scoreEstimation(obs=exceedances, loc=loc, vario=varioModel, weigthFun, dWeigthFun, u = threshold)
  #mvPot::spectralLikelihood(exceedances, loc, varioModel)
}


est <- optim(par =c(1,1),fn = objectiveFunction,
             exceedances = Excedents,loc = LOC,
             vario = f_h_param,weigthFun = weigthFun,
             dWeigthFun = dWeigthFun,threshold =v_seuil,
             control = list(maxit = 4000))
PARAMS_opt<-est$par

variogram_model_score<-function(h){
  #return(abs(h)^1.9)
  return(f_h_param(h,parametre = PARAMS_opt))
}

resultat<-mvPot::simulPareto(n = 1000,loc = LOC,vario =variogram_model_score)
simulations_Pareto<-matrix(unlist(resultat),ncol=37,byrow = TRUE)

matplot(t(simulations_Pareto),type="l")
summary(apply(simulations_Pareto,MARGIN = 1,FUN = calcul_norme_L2))

# cf dombry : multiplier par v seuil --------------------------------------
simulations_Pareto<-simulations_Pareto*v_seuil
matplot(t(simulations_Pareto),type="l")

# Reconversion ------------------------------------------------------------
simulations_Epsi_unif<-(exp(-simulations_Pareto^(-1)))
simulations_EPsi_orig<-simulations_Epsi_unif
for(t in c(1:37)){
  simulations_EPsi_orig[,t]<-((1-simulations_EPsi_orig[,t])^(-DF$forme_t[t])-1)*(DF$échelle_t[t]/DF$forme_t[t])+DF$seuil_t[t]
}
par(mfrow=c(1,2))
matplot(t(simulations_EPsi_orig),type="l")
matplot(t(Residus_donnees[indices_ext_scale_orig,]),type="l")
title("Résultat avec score estimation",outer=TRUE)
par(mfrow=c(1,1))


# Variogram ---------------------------------------------------------------
L<-ncol(Residus_donnees)
vecteur_distances<-c()
Matrice_couples<-list()
z<-1
for(j in 1:L){
  #condition imposée sur le deuxième temps.
  valeurs_t_plus_h<-j:L
  for (i in valeurs_t_plus_h){
    Matrice_couples[[z]]<-c(i,j)
    vecteur_distances<-c(vecteur_distances,abs(i-j))
    z<-z+1
  }
}
Extremogram_donnees<-extremogram_empirique(Matrice_couples=Matrice_couples,
                                           individus_select =Residus_donnees,
                                           Tau=DF$seuil_t)
df<-cbind.data.frame(vecteur_distances,Extremogram_donnees)
colnames(df)<-c("delta","valeur_donnees")
resultat_delta<-df %>% group_by(delta) %>% summarise(val_donnees=mean(valeur_donnees))
formule_variog<-(qnorm((1-(resultat_delta$val_donnees)/2))^2)*2

par(mfrow=c(1,1))
plot(resultat_delta$delta/37,formule_variog, 
     main="Extremogramme empirique")
points(resultat_delta$delta/37,
       variogram_model_score(resultat_delta$delta/37), 
       col="red")


resultat_delta$delta/37

# formule_obtenue ---------------------------------------------------------
estim_modele<-2 * (1 - pnorm(sqrt(variogram_model_score(resultat_delta$delta/37)/ 2)))
plot(resultat_delta$delta,resultat_delta$val_donnees)
points(resultat_delta$delta,estim_modele,col="red")



# Variante ----------------------------------------------------------------
#Compute generating vector
Pareto_obs<-read.csv(file=paste0("residus_clust/",nom_variable,"_obs_ech_Pareto.csv"))[,2:38]
colnames(Pareto_obs)<-c(1:ncol(Pareto_obs))

NL2<-apply(X = Pareto_obs,MARGIN = 1, FUN = calcul_norme_L2)
nu<-quantile(NL2,0.95)
inds_exts<-which(NL2>nu)
Excedents_Censored<-as.list(as.data.frame(t(Pareto_obs[inds_exts,])))


objectiveFunction_Censored = function(parameter, exceedances, loc, vario, threshold){
  if(parameter[1]<0|parameter[2]<0|parameter[1]>2){return(1e50)}
  varioModel <- function(h){
    vario(h = h,param = parameter)
  }
  p <- 499
  latticeRule <- genVecQMC(p, (nrow(LOC) - 1))
  primeP <- latticeRule$primeP
  vec <- latticeRule$genVec
  
  censoredLikelihoodBR(obs =  exceedances, loc = loc, vario =varioModel,
                       u =  threshold, p = primeP, vec = vec)
}


est <- optim(par =c(1,1),fn = objectiveFunction_Censored,
             exceedances = Excedents_Censored,loc = LOC,
             vario = f_h_param,threshold =nu,
             control = list(maxit = 100))
Params_censored<-est$par
Params_censored
variogram_model_score_censored<-function(h){
  #return(abs(h)^1.9)
  return(f_h_param(h,parametre = Params_censored))
}
resultat<-mvPot::simulPareto(n = 1000,loc = LOC,vario =variogram_model_score_censored)
simulations_Pareto<-matrix(unlist(resultat),ncol=37,byrow = TRUE)

matplot(t(simulations_Pareto),type="l")
summary(apply(simulations_Pareto,MARGIN = 1,FUN = calcul_norme_L2))

# cf dombry : multiplier par v seuil --------------------------------------
simulations_Pareto_Censored<-simulations_Pareto*v_seuil
matplot(t(simulations_Pareto_Censored),type="l")

simulations_Epsi_unif<-(exp(-simulations_Pareto_Censored^(-1)))
simulations_EPsi_orig_censored<-simulations_Epsi_unif
for(t in c(1:37)){
  simulations_EPsi_orig_censored[,t]<-((1-simulations_EPsi_orig_censored[,t])^(-DF$forme_t[t])-1)*(DF$échelle_t[t]/DF$forme_t[t])+DF$seuil_t[t]
}
par(mfrow=c(1,2))
matplot(t(simulations_EPsi_orig_censored),type="l")
matplot(t(Residus_donnees[inds_exts,]),type="l")
title("Résultat avec censored",outer=TRUE)
par(mfrow=c(1,1))

plot(resultat_delta$delta,formule_variog)
points(resultat_delta$delta,variogram_model_score(resultat_delta$delta/37),col="green")
points(resultat_delta$delta,variogram_model_score_censored(resultat_delta$delta/37),col="red")

estim_modele_censored<-2 * (1 - pnorm(sqrt(variogram_model_score_censored(resultat_delta$delta/37)/ 2)))
plot(resultat_delta$delta,resultat_delta$val_donnees)
points(resultat_delta$delta,estim_modele_censored,col="red")
points(resultat_delta$delta,estim_modele,col="green")

