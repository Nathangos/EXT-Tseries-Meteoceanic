rm(list=ls())
set.seed(133)
####
M_simulations<-130
#### Changer le nombre 
fnct_p_valeur<-function(M){
  variable_n<-rnorm(M,mean = 0,sd=1)
  a<-cor.test(variable_n,variable_n^2,method = "pearson")
  P1<-a$p.value
  b<-cor.test(variable_n,variable_n^3,method="pearson")
  P2<-b$p.value
  ctest<-cor.test(variable_n^2,variable_n^3,method="pearson")
  P3<-ctest$p.value
  return(c(P1,P2,P3))
}

# On réalise 1000 fois la même expérience ---------------------------------
Resultat<-matrix(replicate(1000,expr = fnct_p_valeur(M=M_simulations)),byrow = TRUE,ncol=3)

# Densité de la p valeur ---------------------------------------
###
plot(density(Resultat[,1]),main="p valeur du test Cov=0 pour C1-C1^2")
abline(v=0.05,col="red")
# Nous avons avec la régression cubique un modèle légèrement significatif. 
###
# Nous pouvons voir si cela peut advenir. 
plot(density(Resultat[,3]),main="p valeur du test Cov=0 pour C1^3-C1^2")
abline(v=0.05,col="red")

# Nombre de rejets de l'hypothèse Cov=0  ------------------------------------------------------
fnct_rejet_pval<-function(vecteur,seuil_p_valeur){
  return(mean(as.numeric(vecteur<seuil_p_valeur)))
}
NB<-apply(X=Resultat,MARGIN = 2,FUN = fnct_rejet_pval,seuil_p_valeur=0.05)
print("probabilitéde rejeter l'hypothèse Cov=0")
NB
# Le résultat pour C1, C1^3 correspond à la réalité -----------------------
# Il peut arriver que C1^2, C1^3 soit corrélé -----------------------------
lambda<-16.69
alpha<-1.32
TEST<-function(x){
  return((x/lambda)^alpha)
}
sqrt(sapply(c(1:37)/37,TEST)*2)

