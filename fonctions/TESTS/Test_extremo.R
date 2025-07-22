source("fonctions/fonctions.R")
require(mev)
alpha<-0.5
beta<-0.8
variog<-function(h){
  return(alpha*(h)^beta)
}
ul<-10
D<-37
vario <- function(x, scale = 0.5, alpha = 0.8){ scale*x^alpha }
Simul_Rpto<-mev::rgparp(n = 1000,
            risk = "l2",
            model = "br",
            vario = vario,
            scale=rep(1,37),
            coord = as.matrix(c(1:37)),
            loc=rep(0,37))
matplot(t(Simul_Rpto),type="l")

Matrice_couples<-list()
L<-ncol(Simul_Rpto)
z<-1
vecteur_distances<-c()
for(j in 1:L){
  #condition imposée sur le deuxième temps.
  valeurs_t_plus_h<-j:L
  for (i in valeurs_t_plus_h){
    Matrice_couples[[z]]<-c(i,j)
    vecteur_distances<-c(vecteur_distances,abs(i-j))
    z<-z+1
  }
}
q_chosen<-0.90
Tau<-sapply(c(1:37),function(x,q,t){
  return(as.numeric(quantile(x[,t],q)))},q=q_chosen,
  x=Simul_Rpto)
Extremogram_simulations<-empirical_extremogram(Matrix_couples =Matrice_couples,Tau = Tau,
                                               inds_select= Simul_Rpto)
Extremogram_simulations
df<-cbind.data.frame(vecteur_distances,Extremogram_simulations)
colnames(df)<-c("delta","data_value")
require(dplyr)
result_delta<-df %>% group_by(delta) %>% 
  summarise(val_data=mean(data_value))
plot(c(0:36),result_delta$val_data)
Result<-sapply(X = c(0:36),
       FUN=function(h){
         2*(1-pnorm((vario(h)/2)^(1/2)))
       })
points(c(0:36),Result,
       col="red")
