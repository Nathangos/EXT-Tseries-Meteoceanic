source("fonctions.R")
f_variogram<-function(s,t){
  coeff<-abs(t-s)
  return(2*(1-exp(-coeff)))
}
f_sigma<-function(t){
  return(2*(1-exp(-abs(t))))
}
set.seed(133)
Realisations_log_norm<-t(Trajectoire_log_norm(f_variogram = f_variogram,f_sigma = f_sigma,Vecteur_temps =c(1:37)/37,M = 1000,alpha = 1))
MHstings<-Procedure_MHastings(Realisations_log_norm,Longueur_echantillon=1000)
matplot(t(Realisations_log_norm),type="l")
Tt<-rbind.data.frame(MHstings)
colnames(Tt)<-1:37
matplot(Tt,type="l")
