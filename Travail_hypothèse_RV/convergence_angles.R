rm(list=ls())
set.seed(133)
par(mfrow=c(1,1))

#### Import des packages. #####
require(dplyr)
require(ggplot2)
require(reshape2)
require(FactoMineR)
require(fda)
require(parallel)
require(mvPot)

### Source des fonctions. ####
source("fonctions_Pareto.R")
source("Modele_ss_tend/Travail_angles.R")
source("fonctions.r")
source("fonctions_r_pareto_EVOLTAU.r")
source("Travail_hypothèse_RV/f_convergence_angle.R")
F_names<-ls()
nb_coeurs<-detectCores()-4


# Ouverture des clusters --------------------------------------------------
coeurs<-makeCluster(nb_coeurs)
clusterExport(coeurs,varlist=F_names)


fonction_evol_quantile<-function(nom_variable,type_donnees,Bornes_Y,N,type_seuils){
    e_<-read.table(paste0("SEUILS_THRESHR_",nom_variable,"_",type_donnees,"_",type_seuils,".txt"))
    Seuils_<-e_$x
    vecteur_k<-seq(50,500,by=1)
    vecteur<-c()
    fonction_export_convergence(coeurs = coeurs,nom = nom_variable,n.dens =N,opt_tend = TRUE,Tau=cummean(Seuils_),BORNES_Y = BORNES_Y,type_donnees = "HIVER")
    OBS_P<-read.csv("OBS_PARETO.csv")[,2:38]
    fnct_k<-function(Obs,k){
      LISTE_p<-fonction_analyse_convergence(Obs=Obs,K= k)
      return(LISTE_p)
    }
    vecteur_CONVERG<-sapply(X = vecteur_k,FUN = fnct_k,Obs=OBS_P)
    MATRICE_moy_coord<-t(cbind.data.frame(vecteur_CONVERG))
    par(mfrow=c(2,2))
    for(j in c(1:ncol(MATRICE_moy_coord))){
      if ((j==1)||(j==5)){
        plot(MATRICE_moy_coord[,j],type="l",ylab=paste0("Moment pour la ",j,"ieme fonction"),xlab="Nombre d'excédents",main=paste0("Premier moment avec ",nom_variable))
        abline(v=110,col="red")
        abline(v=150,col="red")
      }
      else{
        plot(MATRICE_moy_coord[,j],type="l",ylab=paste0("Moment pour la ",j,"ieme fonction"),xlab="Nombre d'excédents")
        abline(v=110,col="red")
        abline(v=150,col="red")
      }
    }
    par(mfrow=c(1,1))
    
    return(MATRICE_moy_coord)
}
for(nom in c("Hs","Surcote","U")){
  Resultat<-fonction_evol_quantile(nom_variable = nom,type_donnees = "HIVER",Bornes_Y = c(0,12),N = 2^(13),type_seuils = "POT")
  
}

# Arret des coeurs --------------------------------------------------------
stopCluster(coeurs)
