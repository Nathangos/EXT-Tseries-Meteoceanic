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
source("Modele_ss_tend/f_ACP_dummy.R")
source("fonctions.r")
source("fonctions_r_pareto_EVOLTAU.r")
F_names<-ls()
nb_coeurs<-detectCores()-4


# Ouverture des clusters --------------------------------------------------
coeurs<-makeCluster(nb_coeurs)
clusterExport(coeurs,varlist=F_names)
## (II) Qualite transformation --------------------------------------------------
#### --------------------------------------------------------------------

dates_Johan<-c("08/03/2008","09/03/2008","10/03/2008")
Param<-c(1,1)

# (I) Resultat des simulations ----------------------------------------------------------------
# #### --------------------------------------------------------------------
Q_H<-0.0008

lancement_code_dummy<-function(nom_variable,Q_H,type_donnees,Bornes_Y){
  
  A<-Simul_dummy(coeurs = coeurs,p_norm = 0.045,nom = nom_variable,opt_tend = TRUE,type_donnees = type_donnees,rho_hissage = Q_H,BORNES_Y = Bornes_Y)
  #Fermeture parallel ------------------------------------------------------
  stopCluster(coeurs)
  return(A)
  
}
R<-lancement_code_dummy(nom_variable = "Hs",Q_H = Q_H,type_donnees = "HIVER",Bornes_Y=c(0,1))
R

# ## ----------------------------------------------------------------------
# Fin du code -------------------------------------------------------------

