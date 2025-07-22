rm(list=ls())

fnct_craft<-function(Vecteur_temps){
  L<-length(Vecteur_temps)
  vecteur_zero<-rep(0,L*L)
  gamma_fonction<-function(s,t){
    h<-abs(t-s)
    return(1-exp(-(h)^(1)))
  }
  f<-SpatialExtremes::covariance(nugget = 0,sill=1,range = 1,smooth = 1,cov.mod = "powexp",plot=FALSE)
  fnct_cov_par_temps<-function(s,t){
    return(f(abs(t-s)))
  }
  Matrice_cov_moitie<-matrix(vecteur_zero,byrow = TRUE,nrow=L)
  for(j in 1:L){
    vecteur_autres_temps<-Vecteur_temps[c(1:j)]
    Ligne<-sapply(vecteur_autres_temps,FUN=fnct_cov_par_temps,s=Vecteur_temps[j])
    Matrice_cov_moitie[j,c(1:j)]<-Ligne    
  }
  Matrice_cov_reste<-t(Matrice_cov_moitie)
  Matrice_cov<-Matrice_cov_moitie+Matrice_cov_reste
  diag(Matrice_cov)<-diag(Matrice_cov)/2
  Realisations_GP<-mvtnorm::rmvnorm(n=10000,sigma = Matrice_cov,method="chol")
  return(Realisations_GP)
}
fnct_SE<-function(Vecteur_temps){
  f_cov<-SpatialExtremes::rgp(n = 10000,coord = Vecteur_temps,nugget = 0,sill=1,range = 1,smooth = 1,cov.mod = "powexp")
  return(f_cov)
}
set.seed(133)
A<-fnct_craft(1:37/37)
par(mfrow=c(2,2))
matplot(t(A),type="l")
A_fun<-apply(A,FUN=var, MARGIN = 2)
A_fun
set.seed(133)
B<-fnct_SE(1:37/37)
matplot(t(B),type="l")
B_fun<-apply(B,FUN=var, MARGIN = 2)
B_fun
