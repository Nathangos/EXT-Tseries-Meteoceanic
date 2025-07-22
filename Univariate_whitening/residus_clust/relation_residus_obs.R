rm(list=ls())
repertoire<-"ss_tend/"
nom_variable<-"Hs"
type_donnees<-"HIVER"
inds_extremes<-paste0("inds_extremes_donnees_",nom_variable,".csv")

source("fonctions.R")
if(type_donnees=="HIVER"){
  repertoire<-paste0(repertoire,"HIVER/")
}
if(type_donnees=="HORS_HIVER"){
  repertoire<-paste0(repertoire,"HORS_HIVER/")
}
# Observer la norme des résidus -------------------------------------------
Residus<-read.csv(paste0("residus_clust/",nom_variable,"_residus.csv"))[,2:38]
colnames(Residus)<-c(1:ncol(Residus))
NL2_residus<-apply(X = Residus,MARGIN = 1,FUN = calcul_norme_L2)

nom_recherche<-paste0(repertoire,nom_variable,"_ss_tend.csv")
repertoire
obs<-read.csv(nom_recherche)[,2:38]
NL2_obs<-apply(X = obs,MARGIN = 1,FUN = calcul_norme_L2)
sub<-which(NL2_residus>0)
plot(NL2_residus[sub],NL2_obs[sub],
     main=paste0("Relation norme epsilon--norme ", nom_variable))
cor.test(NL2_residus,NL2_obs)

# Indices_obs_prises ------------------------------------------------------

fichier_idx_extremes<-paste0("residus_clust/",inds_extremes)
Inds_extremes_L2<-read.csv(fichier_idx_extremes)[,2]
L<-length(NL2_residus)-1
NL2_passe<-NL2_obs[1:L]
NL2_residus_ajd<-NL2_residus[2:length(NL2_residus)]
indice_pos<-which(NL2_residus_ajd>0)

NL2_residus_ajd<-NL2_residus_ajd[indice_pos]
borne_moins<-min(NL2_residus[Inds_extremes_L2])
NL2_passe<-NL2_passe[indice_pos]
print(cor.test(NL2_residus_ajd,NL2_passe))
abline(v=borne_moins,col="red")

# s'intéresse aux individus resid ++ --------------------------------------

indice_part<-which((NL2_residus>5)&(NL2_obs>5))
matplot(t(obs[indice_part,]),type="l")

# Bien une corrélation ----------------------------------------------------

# pour les exts -----------------------------------------------------------
indices_ext<-which(NL2_residus>borne_moins)
cor.test(NL2_residus[indices_ext],NL2_obs[indices_ext])
plot(NL2_residus[indices_ext],NL2_obs[indices_ext], 
     main=paste0("Relation norme epsilon--norme ", nom_variable," pour epsi ext"))


# LM ----------------------------------------------------------------------
modele<-lm(NL2_passe~NL2_residus_ajd)
resume_lm<-summary(modele)

plot(NL2_residus_ajd,NL2_passe, 
     xlab="Norme epsilon", 
     ylab="Norme observation veille")
# Ext et non ext ----------------------------------------------------------
abline(v=borne_moins,col="red")
abline(h=3.5,col="red")
abline(h=2.4,col="orange")
abline(h=3.4,col="orange")
#lines(NL2_residus_ajd,modele$fitted.values,
     # col="green")

modele_non_linear<-nls(formula = NL2_passe~exp(-(Beta0+Beta*NL2_residus_ajd)), 
    start=list(Beta0=0,Beta=1))
objet<-summary(modele_non_linear)
f<-modele_non_linear$m$fitted()
points(NL2_residus_ajd,f,col="green")


# Copules -----------------------------------------------------------------
require(VineCopula)
dataUnif<-pobs(cbind(NL2_residus_ajd,NL2_passe))

Famille_choisie<-BiCopSelect(dataUnif[,1],dataUnif[,2])
summary(Famille_choisie)
Realisations_copules<-BiCopSim(N=1000,family =Famille_choisie$family,par = Famille_choisie$par,
                               par2 =Famille_choisie$par2,check.pars = TRUE)
simulated_resids<-as.numeric(quantile(NL2_residus_ajd,Realisations_copules[,1]))
simulated_passed_NL2<-as.numeric(quantile(NL2_passe,Realisations_copules[,2]))
plot(simulated_resids,simulated_passed_NL2)

modele_cubique<-lm(formula = NL2_passe~I(NL2_residus_ajd^2))
summary(modele_cubique)

elt<-hist(NL2_residus_ajd)

classes<-elt$breaks
l<-length(classes)
sub_cl<-classes[c(2:l)]
ecart_type_categorie<-sapply(X = sub_cl,function(x,categorie,donnees){
  indices<-which((categorie<x)&(categorie>x-0.5))
  return(sd(donnees[indices]))
},categorie=NL2_residus_ajd,donnees=NL2_passe)


plot(elt$mids[c(2:l)],ecart_type_categorie)


# LM+ selection min resid -------------------------------------------------
indices_niv_ext<-which(NL2_residus_ajd>borne_moins)
NL2_resid_ext<-NL2_residus_ajd[indices_niv_ext]
NL2_passe_ext<-NL2_passe[indices_niv_ext]
plot(NL2_resid_ext,NL2_passe_ext)
modele_lm_restricted<-lm(NL2_passe_ext~NL2_resid_ext)
summary(modele_lm_restricted)


# copules -----------------------------------------------------------------
dataUnif2<-pobs(cbind(NL2_resid_ext,NL2_passe_ext))
plot(dataUnif2)
Famille_choisie_tronq<-BiCopSelect(dataUnif2[,1],dataUnif2[,2])
Famille_choisie_tronq

# Independance obtenue ----------------------------------------------------


# Splines -----------------------------------------------------------------
val<-smooth.spline(x = NL2_residus_ajd,y = NL2_passe, 
                   )
plot(NL2_residus_ajd,NL2_passe, 
     xlab="Norme epsilon", 
     ylab="Norme observation veille")
L<-length(NL2_residus_ajd)
ssss<-NL2_passe[2:L]
plot(val$x,ssss)
points(val$x,val$y,col="red",type="l")

summary(val)
help("smooth.spline")

# Choix modele ------------------------------------------------------------
n<-length(NL2_residus_ajd)

# Poids sont constants pour lm et smooth ----------------------------------------------------

RSS_lm<-sum(resume_lm$residuals^2)
print(RSS_lm)
resume_lm$df
BIC_lm<-resume_lm$df[1]*log(n)+n*log(RSS_lm/n)
BIC_lm
RSS_smooth<-val$pen.crit
print(RSS_smooth)
BIC_smooth<-val$df*log(n)+n*log(RSS_smooth/n)
BIC_smooth


