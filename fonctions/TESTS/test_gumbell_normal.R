set.seed(133)
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
source("fonctions.R")
mu<-0
sigma<-2
M<-5000
sortie<-rnorm(M,mean = mu,sd=sigma)
Q1<-as.numeric(quantile(sortie,probs = 0.70))
Q2<-as.numeric(quantile(sortie,probs=0.95))
POT::tcplot(data = sortie,u.range = c(Q1,Q2))
POT::mrlplot(data = sortie,u.range = c(Q1,Q2))
fin<-M/5
vect_k<-c(10:fin)
# Estimateur gamma
####
Graphiques_estimateurs_gamma(series = sortie,vect_k =vect_k ,
                             Titre_graphique = "test",NB_annees = 10)
seuil<-sort(x = sortie,decreasing = TRUE)[250]
seuil
analyse<-fevd(x = sortie,threshold = seuil,type="GP")
resultat<-fonction_ML_extRemes(k = 250,sortie,NB_annees = 10)

