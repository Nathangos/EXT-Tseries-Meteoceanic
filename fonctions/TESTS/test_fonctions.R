rm(list=ls())
set.seed(133)
require(corrplot)
require(dplyr)
source("fonctions.R")

### Test Fréchet. ###
Long_Fr<-2000
shape_<-2
loc_<-1
scale_<-1
vect_Fr<-c(5:2000)
real<-POT::rgpd(Long_Fr,loc=loc_,scale=scale_,shape=shape_)
Vals_hill_pareto<-sapply(vect_Fr,fonction_estimateur_hill,real)
hill(real,"xi")
#avec Frechet. 
real_Frechet<-extRemes::revd(n = Long_Fr,loc = loc_,scale=scale_,shape = shape_,type="GP")
hill(real_Frechet,"xi")
abline(h=shape_)
plot(x=vect_Fr,y=Vals_hill_frech,type="l",xlab="valeur de k",ylab="gamma trouvé",main=paste("Hill plot, alpha=",shape_))

Vals_likelihood_frech<-sapply(vect_Fr,fonction_ML_extRemes,real)
### Attention, le gamma Frechet ne correspond pas exactement au Gamma GEV. ###
fonction_MLplot_resume(resultatML = Vals_likelihood_frech,vecteur_k = vect_Fr,nom_variable = "testFrechet")

# Moments testé. -----------------------------------------------------------
Mtest<-sapply(vect_Fr,fonction_estimateur_moment,real)
plot(Mtest,type="l")

# Frechet -----------------------------------------------------------------

Mtest_frech<-sapply(vect_Fr,fonction_estimateur_moment,real_Frechet)
plot(Mtest_frech,type="l")
# Démonstration : ne pas simuler avec une Fréchet pour tester.
M<-sapply(vect_Fr,POT_fonction_estimateur_moment,real)
plot(M,type="l",col="red")
