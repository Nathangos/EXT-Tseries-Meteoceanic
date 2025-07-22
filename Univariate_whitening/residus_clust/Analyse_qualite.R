rm(list=ls())
set.seed(133)
par(mfrow=c(1,1))

source("../fonctions/fonctions.R")
require(reshape2)
require(ggplot2)
require(gridExtra)
require(grid)
require(dplyr)
require(threshr)
require(extRemes)
require(ggside)
require(kdecopula)

cols_<-c("data"="blue","simulations"="orange","confidence_band"="darkblue")

# Obtenir à partir des simulations  ---------------------------------------
# des résidus de nouvelles trajectoires -----------------------------------
#################

Nom_v<-"Surcote"
Nom_present<-ifelse(Nom_v=="Surcote","Surge",Nom_v)

Temps_interessants<-c(7,13,19,
                      25,31)
inds_extremes<-paste0("residus_clust/inds_extremes_donnees_",Nom_v,".csv")
Modele<-"ACP"
if(Modele=="ACP"){
  simulations<-read.csv(file = paste0("residus_clust/Simulations",Nom_v,"_avec_std.csv"))
}
if(Modele=="BR"){
  simulations<-read.csv(file = paste0("residus_clust/Simulations",Nom_v,"_avec_std_BR.csv"))
}
simulations<-simulations[,c(2:ncol(simulations))]
colnames(simulations)<-c(1:ncol(simulations))


matplot(t(simulations),type="l")
ind_issue<-which.min(simulations[,1])

M<-melt(t(simulations))

colnames(M)<-c("Temps","Individu","valeur")

# Charger le theta -------------------------------------------------------

lien_DF<-paste0("../Travail_hypothèse_RV/EVA_",Nom_v,"_residus.csv")
DF<-read.csv(lien_DF)
print(DF)
M$Theta<-DF$Theta_t[M$Temps]
M$Individu<-as.character(M$Individu)


# Travail sur les simulations vs observations residus ---------------------------------------------
############
Individus_exts<-read.csv(inds_extremes)[,2]

Dates<-read.csv(file="../ss_tend/HIVER/dates_prises.csv")[,2]
d_POIXCT<-as.POSIXct(Dates, format="%d/%m/%Y")
Residus_donnees<-read.csv(file=paste0("residus_clust/",Nom_v,"_residus.csv"))[,c(2:38)]
colnames(Residus_donnees)<-c(1:ncol(Residus_donnees))

# Diagnostic POTS ---------------------------------------------------------
t<-19
Nom_g<-ifelse(Nom_v=="Surcote","S",
              Nom_v)
Outils_POT_graphique(seuil = DF$seuil_t[t],Q1 = 0.50,
                     Q2 = 0.98,series =Residus_donnees[,t],
                     dates = Dates,
                     titre_variable=paste0(Nom_g," at tidal peak (residual)"))



value_kurtosis<-sapply(c(1:ncol(Residus_donnees)),function(x){
  return(mean(((Residus_donnees[,x]-mean(Residus_donnees[,x]))/sd(Residus_donnees[,x]))^(4)))
})
plot(value_kurtosis,main="Evolution aplatissement de la distribution residus")

value_skew<-sapply(c(1:ncol(Residus_donnees)),function(x){
  return(mean(((Residus_donnees[,x]-mean(Residus_donnees[,x]))/sd(Residus_donnees[,x]))^(3)))
})
plot(value_skew,main="Evolution asymetrie de la distribution residus")
# With exponential, it's 6.  ----------------------------------------------


# choice threshold --------------------------------------------------------
t<-19
series_t<-Residus_donnees[,t]
Seuil<-DF$seuil_t[t]
Q1<-0.50
Q2<-0.98
Min_u<-quantile(series_t,Q1)
Max_u<-quantile(series_t,Q2)

POT::mrlplot(data = series_t,u.range = c(Min_u,Max_u),
             main=paste0("Mean residual life plot pour ",Nom_v," en t=",t))
abline(v=Seuil,col="red")

par(mfrow=c(1,2))
POT::tcplot(series_t,ask = FALSE,which = 1,u.range = c(Min_u,Max_u))
abline(v=Seuil,col="red")
POT::tcplot(series_t,ask = FALSE,which = 2,u.range = c(Min_u,Max_u))
abline(v=Seuil,col="red")
par(mfrow=c(1,1))
mtext(paste0("Tcplot pour ",Nom_v," en t=",t),outer=TRUE,line = -2,
      cex=1.3)


Prop_ext<-print(length(Individus_exts)/nrow(Residus_donnees))
Temps<-lubridate::year(d_POIXCT)
Nb_annees<-diff(range(Temps))
NPY<-round(length(Temps)/Nb_annees)

par(mfrow=c(1,2))
matplot(t(Residus_donnees[Individus_exts,]),type="l")
matplot(t(simulations),type="l")
par(mfrow=c(1,1))

M2<-melt(t(Residus_donnees[Individus_exts,]))

colnames(M2)<-c("Temps","Individu","valeur")
M2$Theta<-DF$Theta_t[M2$Temps]
min_y<-min(min(M2$valeur),min(M2$Theta),min(M$valeur))
max_y<-max(max(M2$valeur),max(M2$Theta),max(M$valeur))
M2$Individu<-as.character(M2$Individu)

Titre_EV_simul<-paste0("Simulations/observations obtenues pour la variable ",Nom_v,"(residus) avec limites par (EVA)")
GG1<-ggplot(data=M,aes(x=Temps,y=valeur,group=interaction(Individu),col=Individu))+
  geom_line()+
  ylim(c(min_y,max_y))+
  geom_line(aes(y=Theta,linetype="limite_EV"),col="red")+
  labs(col="Simulation",linetype="Légende",
       caption=paste0("nombre de simulations=",nrow(simulations)))+
  guides(col="none")+
  ggtitle(paste0("Simulations de résidus extrêmes pour ",Nom_v))

GG2<-ggplot(data=M2,aes(x=Temps,y=valeur,group=interaction(Individu),col=Individu))+
  geom_line()+
  geom_line(aes(y=Theta,linetype="limite_EV"),col="red")+
  labs(col="Indice_individu",linetype="Légende")+
  ylim(c(min_y,max_y))+
  ggtitle("Observations extrêmes et limite par EVA")+
  guides(col="none")+
  labs(caption=paste0("nombre d'observations=",length(Individus_exts)))

gridExtra::grid.arrange(GG1,GG2,ncol=2)

colnames(simulations)<-c(1:ncol(simulations))
matplot(t(simulations),type="l")
fnct_quantile<-function(x,q){
  M<-apply(X = x,MARGIN = 2,FUN = function(X){return(as.numeric(quantile(X,q)))})
  return(M)
}

NL2_simulations<-apply(X = simulations,MARGIN = 1,FUN = calcul_norme_L2)

Series_exts<-Residus_donnees[Individus_exts,]
Transpose<-melt(t(Series_exts))
colnames(Transpose)<-c("Time", 
                      "individual", 
                      "value")
Transpose$individual<-as.character(Transpose$individual)
ggplot(data=Transpose,aes(x=Time,y=value,col=individual,
                          group=interaction(individual)))+
  geom_line()+
  xlab("Time (hour) with respect to tidal peak")+
  guides(col="none")
# Percentiles -------------------------------------------------------------

# Simul vs X exts ---------------------------------------------------------
indices_exts_ou_non<-c(1:nrow(Residus_donnees))%in%Individus_exts

vect_percentile<-c(0.05,0.75,0.90,0.95)
resultat<-sapply(vect_percentile,fnct_quantile,x=simulations)
melteur<-melt(t(resultat))
colnames(melteur)<-c("quantite","temps","valeur")
melteur$quantite<-paste0("Q_",vect_percentile[melteur$quantite]*100,"%")
resultat_donnees<-sapply(vect_percentile,fnct_quantile,
                         x=Residus_donnees[Individus_exts,])
N<-nrow(melteur)
melteur_donnees<-melt(t(resultat_donnees))
colnames(melteur_donnees)<-c("quantite","temps","valeur")
melteur_donnees$quantite<-paste0("Q_",vect_percentile[melteur_donnees$quantite]*100,"%")
Percentiles_d_simul<-rbind.data.frame(melteur,melteur_donnees)
Percentiles_d_simul$type<-c(rep("simulations",N),rep("donnees",N))
mu<-colMeans(simulations)

ecart_type<-apply(X = simulations,MARGIN = 2,FUN = sd)
resultat_s<-melt(t(cbind.data.frame(mu,ecart_type)))
colnames(resultat_s)<-c("quantite","temps","valeur")
mu<-colMeans(Residus_donnees[Individus_exts,])
ecart_type<-apply(X =Residus_donnees[Individus_exts,],MARGIN = 2,FUN = sd)
resultat_d<-melt(t(cbind.data.frame(mu,ecart_type)))
colnames(resultat_d)<-c("quantite","temps","valeur")

ggplot(data=Percentiles_d_simul,aes(x=temps,y=valeur,col=type))+
  geom_line()+
  facet_wrap(~quantite)+
  ggtitle(paste0("Percentiles données vs simulations pour ",Nom_v," (résidus)"))+
  labs(caption = paste0("nombre_simulations=",nrow(simulations),", nombre d'individus=",length(Individus_exts)))

R_mu_sd<-rbind(resultat_s,resultat_d)
R_mu_sd$type<-c(rep("simulations",nrow(resultat_d)),rep("donnees",nrow(resultat_d)))
ggplot(data=R_mu_sd,aes(x=temps,y=valeur,col=type))+
  geom_line()+
  facet_wrap(~quantite,scales = "free_y")+
  ggtitle(paste0("Moments données vs simulations pour ",Nom_v," (résidus)"))+
  labs(caption = paste0("nombre_simulations=",nrow(simulations),", nombre d'individus=",length(Individus_exts)))

# Extremogramme -----------------------------------------------------------
L<-ncol(Residus_donnees[Individus_exts,])
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


Tau<-sapply(c(1:37),function(x,q,t){
  return(as.numeric(quantile(x[,t],q)))},q=0.90,x=Residus_donnees[Individus_exts,])
Tau2<-sapply(c(1:37),function(x,q,t){
  return(as.numeric(quantile(x[,t],q)))},q=0.90,x=simulations)
Extremogram_donnees<-empirical_extremogram(Matrix_couples=Matrice_couples,
                                          inds_select  =Residus_donnees[Individus_exts,],
                                           Tau = Tau)
Extremogram_simul<-empirical_extremogram(Matrix_couples=Matrice_couples,
                                         inds_select  =simulations,
                                         Tau = Tau)
type_corr<-"pearson"
Coor_lien_real<-estim_Pearson_tidal_cycle(Matrice_couples = Matrice_couples,
                          individus_select = Residus_donnees[Individus_exts,],
                          type_corr = "pearson")
Coor_lien_simul<-estim_Pearson_tidal_cycle(Matrice_couples = Matrice_couples,
                                          individus_select = simulations,
                                          type_corr ="pearson")

df_coor_lien<-cbind.data.frame(Coor_lien_real,Coor_lien_simul,
                               vecteur_distances)
colnames(df_coor_lien)<-c("valeur_simulation","valeur_donnees",
                          "delta")
resultat_corr<-df_coor_lien %>% group_by(delta) %>% summarise(val_donnees=mean(valeur_donnees),
                                                     val_simul=mean(valeur_simulation))
ggplot(data=resultat_corr,aes(x=delta,y=val_simul,col="simulations"))+
  geom_point()+
  geom_point(aes(y=val_donnees,col="data"))+
  scale_color_manual(values = cols_)+
  ggtitle(paste0("Simul vs donnees pour la corr ", type_corr))+
  ylab("value")+xlab("Lag_h")

df<-cbind.data.frame(Extremogram_simul,Extremogram_donnees,vecteur_distances)
colnames(df)<-c("valeur_simulation","valeur_donnees","delta")

resultat_delta<-df %>% group_by(delta) %>% summarise(val_donnees=mean(valeur_donnees),
                                                     val_simul=mean(valeur_simulation))

Resultat<-replicate(n = 500,expr = fnct_estim_extremo_resample(B = 130,inds_ext  = Residus_donnees[Individus_exts,],
                                                   Tau = Tau, 
                                                   Matrix_couples  = Matrice_couples, 
                                                   vector_distances = vecteur_distances))

Q05<-fnct_quantile(x = t(Resultat),q=0.05)
Q95<-fnct_quantile(x = t(Resultat),q=0.95)
resultat_delta$borne_inf<-Q05
resultat_delta$borne_sup<-Q95
Q05
travail_extremo<-as.data.frame(resultat_delta[,c("val_donnees","val_simul","borne_inf","borne_sup")])
travail_extremo$Temps<-c(1:nrow(travail_extremo))

GG_extremo<-ggplot(travail_extremo,aes(x=Temps,y=val_donnees,col="data"))+
  geom_line()+
  geom_point()+
  geom_line(aes(y=val_simul,col="simulations"))+
  geom_point(aes(y=val_simul,col="simulations"))+
  geom_ribbon(mapping = aes(ymin=borne_inf,ymax=borne_sup,col="confidence_band"),alpha=0.15,
              fill="grey", linetype = "dashed")+
  scale_color_manual(values=cols_)+
  ggtitle(paste0("Estimated extremogram for ",Nom_v," (residual) with resampling"))+
  labs(caption=paste0("Nb_samples=",500,", Size sample=",130,
                      ", quantile_q=",0.90),
       col="Légende")

print(GG_extremo)

# Compar_results_BR -------------------------------------------------------
BR_test<-read.csv(file="residus_clust/SimulationsSurcote_avec_std_BR.csv")
simulations_BR<-BR_test[,c(2:38)]
colnames(simulations_BR)<-c(1:37)
matplot(t(simulations_BR),type="l")
Tau<-sapply(c(1:37),function(x,q,t){
  return(as.numeric(quantile(x[,t],q)))},q=0.90,x=Residus_donnees[Individus_exts,])
Tau2<-sapply(c(1:37),function(x,q,t){
  return(as.numeric(quantile(x[,t],q)))},q=0.90,x=simulations_BR)
Extremogram_donnees<-empirical_extremogram(Matrix_couples=Matrice_couples,
                                           inds_select  =Residus_donnees[Individus_exts,],
                                           Tau = Tau)
Extremogram_simul_BR<-empirical_extremogram(Matrix_couples=Matrice_couples,
                                         inds_select  =simulations_BR,
                                         Tau = Tau)
df_BR<-cbind.data.frame(Extremogram_simul_BR,Extremogram_donnees,vecteur_distances)
colnames(df_BR)<-c("valeur_simulation_BR","valeur_donnees","delta")

resultat_delta_BR<-df_BR %>% group_by(delta) %>% summarise(val_donnees=mean(valeur_donnees),
                                                     val_simul=mean(valeur_simulation_BR))
resultat_delta_BR$borne_inf<-Q05
resultat_delta_BR$borne_sup<-Q95
travail_extremo_BR<-travail_extremo
travail_extremo_BR<-as.data.frame(resultat_delta_BR[,c("val_donnees","val_simul","borne_inf","borne_sup")])
travail_extremo_BR$Temps<-c(1:nrow(travail_extremo))
GG_extremo<-ggplot(travail_extremo_BR,aes(x=Temps,y=val_donnees,col="data"))+
  geom_line()+
  geom_point()+
  geom_line(aes(y=val_simul,col="simulations"))+
  geom_point(aes(y=val_simul,col="simulations"))+
  geom_ribbon(mapping = aes(ymin=borne_inf,ymax=borne_sup,col="confidence_band"),alpha=0.15,
              fill="grey", linetype = "dashed")+
  scale_color_manual(values=cols_)+
  ggtitle(paste0("Estimated extremogram (BR) for ",Nom_v," (residual) with resampling"))+
  labs(caption=paste0("Nb_samples=",500,", Size sample=",130,
                      ", quantile_q=",0.90),
       col="Légende")

print(GG_extremo)

# Analyse norme L2 --------------------------------------------------------

NL2_observations<-apply(X = Residus_donnees,MARGIN = 1,FUN = calcul_norme_L2)

df_norme<-cbind.data.frame(Temps,NL2_observations)
colnames(df_norme)<-c("time","obs")
par(mfrow=c(2,2))
Lims<-c(quantile(df_norme$obs,0.50),
        quantile(df_norme$obs,0.98))
POT::tcplot(data = NL2_observations,u.range = Lims)
POT::mrlplot(data = NL2_observations,u.range=Lims)
POT::diplot(data = df_norme,u.range = Lims,nt=300,main = paste0("Analyse diplot norme de ",Nom_v))
abline(h=1,col="red")
mtext("test",outer=TRUE,line=-1)
par(mfrow=c(1,1))


# # Travail sur la norme ----------------------------------------------------
# #############################

periods_years<-c(2,5,10,20,50,80,
                 100,120,200,
                 250,300,500,
                 800)
NB_annee<-37
limite_sup<-3*37
periods_years<-periods_years[which(periods_years<limite_sup)]

# 1°) loi avec orig et avec négatif ---------------------------------------
seuil<-quantile(x = NL2_observations,0.95)
NPY<-round(length(NL2_observations)/37)
modele_ev<-extRemes::fevd(x = NL2_observations,threshold = seuil,
                          time.units = paste0(NPY,"/year"),
                          type="GP")
modele_ev$results$par
RL_ggplot(series = NL2_observations,seuil = seuil,
          period_years = periods_years,
          NPY = NPY,titre = paste0("EVA of the L2 norm of ",Nom_v),
          nom_variable = Nom_v, 
          unit_used = "m")

# 2° variante, enlever si obs negative ------------------------------------
Indicatrice_neg_un_temps<-apply(Residus_donnees<0,MARGIN = 1,FUN = sum)
Obs_non_neg<-Residus_donnees[which(Indicatrice_neg_un_temps==0),]
norme_diff<-apply(Obs_non_neg, MARGIN=1, FUN=calcul_norme_L2)
seuil_diff<-quantile(norme_diff,0.90)
RL_ggplot(series =norme_diff,seuil = seuil_diff,
          period_years = periods_years,
          unit_used = "m",
          NPY = NPY,titre = "test",nom_variable = Nom_v)

NL2_exts<-NL2_observations[Individus_exts]
indice_min<-which(NL2_exts==min(NL2_exts))
indice_base_orig<-Individus_exts[indice_min]
y<-Residus_donnees[indice_base_orig,]
# EVA-donnees -------------------------------------------------------------

# Mettre à jour le seuil utilisé-------
# Faire tourner à mon avis ----------------
# Problème vient du passage à l'échelle orig ----------------------------

seuil<-min(NL2_exts)

indicatrice_ext_ou_non<-c(1:length(NL2_observations))%in%Individus_exts
indicatrice_ext_ou_non<-ifelse(indicatrice_ext_ou_non,"ext","non_ext")
df_ind_niveau<-cbind.data.frame(NL2_observations,indicatrice_ext_ou_non)
colnames(df_ind_niveau)<-c("val_L2","indicatrice")
ggplot(data=df_ind_niveau,aes(y=val_L2,x=indicatrice,fill=indicatrice))+
  geom_violin()+
  geom_boxplot(width=0.2)+
  ylab("L2 norm")+
  geom_hline(aes(color="seuil_L2",yintercept = seuil))+
  labs(col="Légende")+
  ggtitle(paste0("Densité de la norme L2 selon le type d'individu ",Nom_v))

# proportion_donnees -------------------------------------------------------
inds_exts_pbs<-which((NL2_observations>seuil)&(indices_exts_ou_non==FALSE))

matplot(t(Residus_donnees[inds_exts_pbs,]),type="l")
Mat<-Residus_donnees[inds_exts_pbs,]>0
resultat<-as.numeric(apply(X = Mat,MARGIN = 1,FUN = sum))
positive<-which(resultat==ncol(Residus_donnees))

fusion<-c(inds_exts_pbs[positive],indice_base_orig)

ss_base_pb<-Residus_donnees[fusion,]
matrice<-melt(t(ss_base_pb))
colnames(matrice)<-c("time","index","value")
matrice$ds_base_ext<-ifelse(matrice$index==indice_base_orig,"extreme","non_extreme")
matrice$index<-as.character(matrice$index)
ggplot(data=matrice,aes(x=time,y=value,linetype=ds_base_ext,group=interaction(index),
                        col=index))+
  geom_line()+
  xlab("Time (hour) with respect to tidal peak")+
  guides(col="none")+
  ggtitle(paste0("Présentation du cas particulier"))
  
mtext("Matplot des observations ds la queue avec min et celles proches (echelle originale)")

Pb<-NL2_observations[inds_exts_pbs[positive]]
summary(Pb)

# on enleve les traject negatives -----------------------------------------

col<-Dates[inds_exts_pbs]
valeur<-NL2_observations[inds_exts_pbs]
data_illustr<-cbind.data.frame(valeur,col,c(1:length(valeur)))
colnames(data_illustr)<-c("valeur","date","indice")

cas_particuliers<-Residus_donnees[inds_exts_pbs,]

modele_ML<-fevd(x = NL2_observations,threshold = seuil,type="GP",
                time.units = paste0(NPY,"/year"))
NL2_au_dessus<-subset(NL2_observations,NL2_observations>seuil)
plot(modele_ML)
gamma_obs<-modele_ML$results$par[[2]]
sigma_obs<-modele_ML$results$par[[1]]

print("Estimateurs ML pour les observations")
print(c(gamma_obs,sigma_obs))

AD<-goftest::ad.test(x=NL2_simulations,null=extRemes::"pevd",
                 threshold=seuil,
                 shape=gamma_obs,
                 scale=sigma_obs,type="GP")
KS<-ks.test(x=NL2_simulations,extRemes::"pevd",
        threshold=seuil,
        shape=gamma_obs,
        scale=sigma_obs,type="GP")
print(c(AD$p.value,KS$p.value))
plot(density(NL2_simulations))
plot(density(NL2_au_dessus))


# essai avec MAJ des parametres ------------------------------------------


# EVA-simulations ---------------------------------------------------------------------

# estimateur<-fonction_ML_extRemes(k = 800,donnees = NL2_simulations,
#                      typeML = "GP",NB_annees = 37 )
# Graphiques_estimateurs_gamma(series = NL2_simulations,vect_k = c(10:600),
#                              Titre_graphique = paste0("Estimator of the shape parameter of ",Nom_v,"' norm"),NB_annees = 37)
# 
npy_orig<-round(length(NL2_observations)/37)



# Pareto ------------------------------------------------------------------

Pareto_obs<-read.csv(file=paste0("residus_clust/Pareto/",Nom_v,"_obs_ech_Pareto.csv"))[,c(2:38)]
colnames(Pareto_obs)<-c(1:ncol(Pareto_obs))

# L2_Pareto ---------------------------------------------------------------
L2_PTO<-apply(X = Pareto_obs,MARGIN = 1,FUN = calcul_norme_L2)
Mu<-35
Resultat_threshr<-Outils_Threshr_choix(Q1 = 0.90,
                     Q2 = 0.98,
                     NT_ths = 50,
                     variable = L2_PTO,
                     N_v = 2)
plot(Resultat_threshr)

Outils_POT_graphique(seuil =Mu,
                     Q1 = 0.90,
                     Q2 = 0.98,
                     series = L2_PTO,
                     dates =Dates,
                     titre_variable = "Norm for tranf")

# Evol gamma with other methods -------------------------------------------

vect_k<-c(20:600)
Nom_graph<-ifelse(Nom_v=="Surcote","S",Nom_v)
TITLE_graph<-paste0("Shape parameter T(",Nom_graph,")'s L2 norm (residual)")

L2_orig<-apply(X = Residus_donnees,MARGIN = 1,FUN = calcul_norme_L2)
GG_ext1<-Graphics_estimators_gamma(series = L2_orig,
                          vect_k = vect_k,
                          Title_graphic = NULL,
                          NB_years = NULL)

GG_ext2<-Graphics_estimators_gamma(series = L2_PTO,
                             vect_k = vect_k,
                          Title_graphic = NULL,
                             NB_years = NULL)
All_exts<-GG_ext1+GG_ext2
All_exts
ggsave(filename = "graphiques_Surcote/Shape_parameter_S_TS.png",
       plot = All_exts,
       width = 12,
       height = 5)

Tau<-apply(X = Pareto_obs,MARGIN = 2,FUN = function(x){quantile(x,0.95)})
Extremo_<-empirical_extremogram(Matrix_couples=Matrice_couples,
                      inds_select =Pareto_obs[Individus_exts,],
                      Tau = Tau)
vecteur_distances
df<-cbind.data.frame(Extremo_,vecteur_distances)
colnames(df)<-c("val_donnees","delta")

SS_ech_gnal<-sample(c(1:nrow(Residus_donnees)),size = 20)
Echs<-melt(t(Residus_donnees[SS_ech_gnal,]))
colnames(Echs)<-c("Time", "individual", 
                         "value")
Q95<-apply(Residus_donnees,MARGIN = 2,FUN = function(x){return(quantile(x,0.95))})
Q05<-apply(Residus_donnees,MARGIN = 2,FUN = function(x){return(quantile(x,0.05))})
Mean<-colMeans(Residus_donnees)
Mean
Echs$individual<-as.character(Echs$individual)
Echs$Q05<-Q05[Echs$Time]
Echs$Q95<-Q95[Echs$Time]

# Chargement_EV -----------------------------------------------------------

resultat<-read.csv(file=paste0("Travail_hypothèse_RV/EVA_",Nom_v,"_residus.csv"))
mu<-resultat$seuil_t
Echs$seuil<-mu

#titre_ech<-expression("Sampling of 20 observations of "~epsilon[M]~" for "~Nom_v)

ggplot(Echs,aes(x=Time,y=value,col=individual,
                       group=interaction(individual)))+
  geom_ribbon(mapping=aes(ymin=Q05,ymax=Q95,linetype="confidence_band"),
              alpha=0.02,
              fill="grey",col="darkblue")+
  ylab("value (m)")+
  guides(col="none")+
  xlab("Time (hour) with respect to tidal peak")+
  geom_line(aes(y=seuil,linetype="seuil"),col="red")+
  geom_line(alpha=0.7)+
  theme(axis.title=element_text(size=15))+
  scale_linetype_manual("Legend",values=c("confidence_band"=2, 
                                          "seuil"=5))
  #ggtitle(label = do.call("substitute", list(titre_ech[[1]], 
  #                                           list(Nom_v = Nom_present))))



Pareto_simul<-read.csv(file=paste0("residus_clust/Pareto/",Nom_v,"_simul_Pareto.csv"))
Pareto_simul_BR<-read.csv(file=paste0("residus_clust/Pareto/",Nom_v,"_simul_BR_Pareto.csv"))

Pareto_simul<-Pareto_simul[,c(2:38)]
Resultat_FONCTION<-apply(Pareto_simul,MARGIN = 2,FUN = function(x){
  return(extRemes::pevd(q = x,threshold = 1,scale = 1,
                        shape = 1,type = "GP"))
})
apply(Resultat_FONCTION,MARGIN = 2,FUN = min)


Pareto_angle<-t(t(Pareto_obs)%*%diag(apply(MARGIN = 1,
                                         FUN = calcul_norme_L2,
                                         X=Pareto_obs)^(-1)))
min(Pareto_angle[,1])

ss_ech<-sample(Individus_exts,size = 20,replace = FALSE)
Q95_orig<-apply(Residus_donnees[Individus_exts,],MARGIN = 2,FUN = function(x){return(quantile(x,0.95))})
Q05_orig<-apply(Residus_donnees[Individus_exts,],MARGIN = 2,FUN = function(x){return(quantile(x,0.05))})
Mean_orig<-colMeans(Residus_donnees[Individus_exts,])

repres_Orig<-melt(t(Residus_donnees[ss_ech,]))
colnames(repres_Orig)<-c("Time", 
                           "individual", 
                           "value")
repres_Orig$Q95<-Q95_orig[repres_Orig$Time]
repres_Orig$mean<-Mean_orig[repres_Orig$Time]
repres_Orig$Q05<-Q05_orig[repres_Orig$Time]
repres_Orig$individual<-as.character(repres_Orig$individual)
repres_Orig$Time<-(repres_Orig$Time-19)/6
titre_orig<-expression("Sampling of 20 extremes observations of "~epsilon[M]~" for "~Nom_v)

GG_orig<-ggplot(repres_Orig,aes(x=Time,y=value,col=individual,
                         group=interaction(individual)))+
  geom_ribbon(mapping=aes(ymin=Q05,ymax=Q95,linetype="confidence_band"),
              alpha=0.02,
              fill="grey",col="darkblue")+
  ylab(paste0(Nom_present, " (m)"))+
  geom_line(aes(y=mean,linetype="mean"),col="red")+
  geom_line(alpha=0.7)+
  guides(col="none")+
  theme(axis.title=element_text(size=15),
        legend.text=element_text(size=10))+
  theme(axis.title.x = element_blank())+
  scale_linetype_manual("Legend",values=c("confidence_band"=2,"mean"=4))
  #ggtitle(label = do.call("substitute", list(titre_orig[[1]], 
  #                                           list(Nom_v = Nom_present))))

repres_Pareto<-melt(t(Pareto_obs[ss_ech,]))

colnames(repres_Pareto)<-c("Time", 
                       "individual", 
                       "value")
repres_Pareto$individual<-as.character(repres_Pareto$individual)
titre<-expression("Sampling of 20 extremes observations of T("~epsilon[M]~") for "~Nom_v)
Nom_present<-ifelse(Nom_v=="Surcote","Surge",Nom_v)
Q95<-apply(Pareto_obs[Individus_exts,],MARGIN = 2,FUN = function(x){return(quantile(x,0.95))})
Q05<-apply(Pareto_obs[Individus_exts,],MARGIN = 2,FUN = function(x){return(quantile(x,0.05))})
Mean_Pareto<-colMeans(Pareto_obs[Individus_exts,])
repres_Pareto$Q95<-Q95[repres_Pareto$Time]  
repres_Pareto$Q05<-Q05[repres_Pareto$Time]
repres_Pareto$mean<-Mean_Pareto[repres_Pareto$Time]
repres_Pareto$Time<-(repres_Pareto$Time-19)/6
GG_Pareto<-ggplot(repres_Pareto,aes(x=Time,y=value,col=individual,
                          group=interaction(individual)))+
  geom_ribbon(mapping=aes(ymin=Q05,ymax=Q95,linetype="confidence_band"),alpha=0.02,
              fill="grey",col="darkblue")+
  geom_line(aes(y=mean,linetype="mean"),col="red")+
  geom_line(alpha=0.7)+
  guides(col="none")+
  ylab(paste0("T(",Nom_present,") (-)"))+
  theme(axis.title=element_text(size=15),
        legend.text=element_text(size=10))+
  theme(axis.title.x = element_blank())+
  scale_linetype_manual("Legend",values=c("confidence_band"=2,"mean"=4))
  #ggtitle(label = do.call("substitute", list(titre[[1]], 
  #                                         list(Nom_v = Nom_present))))

shared_xlab <- grid::textGrob("Time (hour) with respect to tidal peak", gp = grid::gpar(fontsize = 15))
combined<-(GG_orig|GG_Pareto) 
require(patchwork)
All_orig_ext<-patchwork::wrap_elements(full = combined) / patchwork::wrap_elements(shared_xlab) +
  plot_layout(heights = c(1, 0.05))
ggsave("graphiques_Surcote/exs_S.png",
       plot = All_orig_ext,
       width=12,
       height=7)

# Xbottom<-textGrob(,
#                   gp = gpar(col = "black", fontsize = 15))
# grid.arrange(GG_orig,GG_Pareto,
#              bottom=Xbottom,ncol=2)
L2_Pareto_obs<-apply(X = Pareto_obs,MARGIN = 1,FUN=calcul_norme_L2)

# Travail sur la simulation de nouveaux temps d'exts ----------------------
# -------------------------------------------------------------------------


# # Loi du maximum ds_donnees -----------------------------------------------
number_years<-37
npy<-round(length(L2_Pareto_obs)/number_years)
npy



# Travail sur elt ---------------------------------------------------------
obs_exts<-Residus_donnees[Individus_exts,]
DF

#Temps_interessants

for(t in c(19)){
  quantile_t<-1-DF$p_u_t[t]
  series_simul<-simulations[,t]
  series_t<-Residus_donnees[,t]

  # Opitz_ML ----------------------------------------------------------------
  l_t<-list("shape"=DF$forme_t[t],"scale"=DF$échelle_t[t],
            "threshr"=DF$seuil_t[t])
  seuil<-quantile(series_t,quantile_t)
  lag<-t-19
  print(c(min(subset(series_simul,series_simul>seuil)),seuil))
  comment<-ifelse(lag<0,paste0(abs(round(lag/6)),"h before"),ifelse(lag==0," at",paste0(round(lag/6),"h after")))
  
  # resultat<-RL_ggplot(series = series_t,seuil = seuil,
  #           period_years = periods_years,
  #           NPY = npy_orig,
  #           plus_simul = TRUE,
  #           series_simul = series_simul, 
  #           titre = paste0("Return level (ML) ",comment, " the tidal peak for ",Nom_v," residual"), 
  #           nom_variable = Nom_v, 
  #           cols_ggplot = cols_, 
  #           unit_used = "m",
  #           alpha=0.01)
  # 
  # resultat_2opitz<-RL_ggplot_boot_opitz(series = series_t,seuil = seuil,
  #                                       period_years = periods_years,
  #                                       NPY = npy_orig,
  #                                       titre = paste0("Return level (Opitz_ML) ",comment, " the tidal peak for ",Nom_v," residual"),
  #                                       nom_variable = Nom_v,
  #                                       liste_ML =  l_t,
  #                                       B =250,series_simul = series_simul,
  #                                       alpha = 0.01,unit_used="m",
  #                                       cols_ggplot = cols_)
}


# logique à l'oeuvre : condition l intervient ds le qqplot ----------------
require(scales)
resultat<-RL_ggplot_cond_ext(series = series_t,seuil = seuil,
          period_years = periods_years,
          NPY = npy_orig,
          plus_simul = TRUE,
          series_simul = series_simul, 
          titre = paste0("Return level adjusted (ML) ",comment, " the tidal peak for ",Nom_v," residual"), 
          nom_variable = Nom_v, 
          cols_ggplot = cols_, 
          unit_used = "m",
          alpha=0.01,Individus_exts = Individus_exts)

# Forme des series temporelles --------------------------------------------
# Work on independence ----------------------------------------------------


# original scale ----------------------------------------------------------
NL2<-apply(X = obs_exts,MARGIN = 1,FUN = calcul_norme_L2)
Theta<-t(t(obs_exts)%*%diag(NL2^(-1)))
matplot(t(Theta),type="l")
Analyse_theta<-FactoMineR::PCA(X = Theta,scale.unit =TRUE,graph = FALSE)
lambdos<-Analyse_theta$eig[,1]
Coords<-Analyse_theta$ind$coord[,1:3]

plot(1-c(0,cumsum(lambdos)/sum(lambdos)))

# elbow criteria-----------------------------------------------------------
Coords_acp<-Analyse_theta$ind$coord[,c(1:3)]
liste_indep_Angle_orig<-matrix(NA,nrow=ncol(Coords_acp),
                                    ncol=3)
for(j in c(1:ncol(Coords_acp))){
  indice<-j
  liste_indep_Angle_orig[j,]<-
    c(cor.test(NL2,Coords_acp[,indice],
               method="pearson")$p.value,
      cor.test(NL2,Coords_acp[,indice],
               method="spearman")$p.value,
      cor.test(NL2,Coords_acp[,indice],
               method="kendall")$p.value)
}
colnames(liste_indep_Angle_orig)<-c("pearson",
                                         "spearman",
                                         "kendall")
liste_indep_Angle_orig

# ech pareto --------------------------------------------------------------
coords_donnees_Pareto<-read.csv(file=paste0("residus_clust/ACP_results/",Nom_v,"_coords_obs.csv"))
coords_donnees_Pareto<-coords_donnees_Pareto[,-1]
L2_norme_exts<-apply(Pareto_obs[Individus_exts,],
                     MARGIN = 1,FUN = calcul_norme_L2)
liste_indep_Angle_Pareto<-matrix(NA,nrow=ncol(coords_donnees_Pareto),
                            ncol=3)
for(j in c(1:ncol(coords_donnees_Pareto))){
  indice<-j
  liste_indep_Angle_Pareto[j,]<-
    c(cor.test(L2_norme_exts,coords_donnees_Pareto[,indice],
                                   method="pearson")$p.value,
      cor.test(L2_norme_exts,coords_donnees_Pareto[,indice],
               method="spearman")$p.value,
      cor.test(L2_norme_exts,coords_donnees_Pareto[,indice],
               method="kendall")$p.value)
}
colnames(liste_indep_Angle_Pareto)<-c("pearson",
                                         "spearman",
                                         "kendall")

liste_indep_Angle_Pareto
liste_indep_Angle_orig

# Analyse scores ----------------------------------------------------------
nom_present<-ifelse(Nom_v=="Surcote","S",Nom_v)
nom_fichier<-paste0("residus_clust/ACP_results/Theta_coords_ACP_",Nom_v,".csv")
nom_fichier
nom_fi_simul<-paste0("residus_clust/ACP_results/Theta_Vine_simul_",Nom_v,".csv")
resultat_obs<-read.csv(nom_fichier)[,-1]

resultat_simul<-read.csv(nom_fi_simul)[,-1]
colnames(resultat_obs)<-paste0("C_",seq.int(from = 1,to = ncol(resultat_obs)))
colnames(resultat_simul)<-colnames(resultat_obs)
ks.test(resultat_obs[,1],resultat_simul[,1])
ks.test(resultat_obs[,2],resultat_simul[,2])
ks.test(resultat_obs[,3],resultat_simul[,3])

# Import inertia ----------------------------------------------------------
Table_inertia<-read.csv(paste0("residus_clust/ACP_results/evolution_eigen_value_",Nom_v,".csv"))
Prop_inertia<-round(Table_inertia[,2]/sum(Table_inertia[,2]),2)*100
resultat_sample<-resultat_simul

Titre_theta<-expression("Scores on a PCA basis of "~Theta~" for "~Nom_present)

df_combo<-rbind.data.frame(resultat_obs,resultat_sample)
df_combo$Legend<-c(rep("data",nrow(resultat_obs)),
                 rep("simulations",nrow(resultat_sample)))

expr1<-function_expression_prop_variance_j(j = 1,prop_variance = Prop_inertia)
expr2<-function_expression_prop_variance_j(j = 2,prop_variance = Prop_inertia)

expr1
GG1<-ggplot(data=df_combo,aes(x=C_1,y=C_2,colour=Legend,shape=Legend),
)+
  geom_point(aes(size=Legend))+
  geom_xsidedensity(data=df_combo,aes(fill=Legend), alpha = 0.5)+
  geom_ysidedensity(data=df_combo,aes(fill=Legend), alpha = 0.5)+
  xlab(expr1)+
  ylab(expr2)+
  scale_color_manual(values=cols_)+
  scale_fill_manual(values=cols_)+
  scale_shape_manual(values = c("simulations"=17,"data"=19))+
  scale_size_manual(values=c("simulations"=0.75,"data"=1.5))+
  guides(fill="none")+
  theme(axis.title=element_text(size=15),
        legend.text=element_text(size=10))

GG1

if(ncol(resultat_sample)>2){
  # (C_2,C_3)
  expr_<-expression(C[3]~" ("~prop3~"% of the variance)")
  expr_<-do.call("substitute", 
                 list(expr_[[1]], list(prop3 = Prop_inertia[3])))

  GG2<-ggplot(data=df_combo,aes(x=C_2,y=C_3,colour=Legend,shape=Legend),
  )+
    geom_point(aes(size=Legend))+
    geom_xsidedensity(data=df_combo,aes(fill=Legend), alpha = 0.5)+
    geom_ysidedensity(data=df_combo,aes(fill=Legend), alpha = 0.5)+
    xlab(expr2)+
    ylab(expr_)+
    scale_color_manual(values=cols_)+
    scale_fill_manual(values=cols_)+
    scale_shape_manual(values = c("simulations"=17,"data"=19))+
    scale_size_manual(values=c("simulations"=0.75,"data"=1.5))+
    guides(fill="none")+
    theme(axis.title=element_text(size=15),
          legend.text=element_text(size=10))

  
  # (C_1,C_3)
  GG3<-ggplot(data=df_combo,aes(x=C_1,y=C_3,colour=Legend,shape=Legend),
  )+
    geom_point(aes(size=Legend))+
    geom_xsidedensity(data=df_combo,aes(fill=Legend), alpha = 0.5)+
    geom_ysidedensity(data=df_combo,aes(fill=Legend), alpha = 0.5)+
    xlab(expr1)+
    ylab(expr_)+
    scale_color_manual(values=cols_)+
    scale_fill_manual(values=cols_)+
    scale_shape_manual(values = c("simulations"=17,"data"=19))+
    scale_size_manual(values=c("simulations"=0.75,"data"=1.5))+
    guides(fill="none")+
    theme(axis.title=element_text(size=15),
          legend.text=element_text(size=10))
  # GG3<-ggplot(data=df_combo,aes(x=C_1,y=C_3,colour=Legend,shape=Legend,
  #                               size=Legend),
  # )+
  #   geom_point(aes(size=Legend))+
  #   geom_xsidedensity(data=df_combo,aes(fill=Legend), alpha = 0.5)+
  #   geom_ysidedensity(data=df_combo,aes(fill=Legend), alpha = 0.5)+
  #   xlab(expr1)+
  #   ylab(expr_)+
  #   #labs(colour="Legend")+
  #   scale_color_manual(values=cols_)+
  #   scale_fill_manual(values=cols_)+
  #   scale_shape_manual(values = c("simulations"=17,"data"=19))+
  #   scale_size_manual(values=c("simulations"=0.75,"data"=1.5))+
  #   guides(fill="none")+
  #   theme(axis.title=element_text(size=15))
  # #ggtitle(label = do.call("substitute", list(Titre_theta[[1]], 
  # #                                           list(Nom_present = Nom_present))))
  
}
GG2
GG3

blank <- ggplot(width=0.5) + theme_void()

# layout <- c(
#   area(1, 1, 1, 1),  # GG1 in row 1, col 1
#   area(1, 2, 1, 2),  # GG2 in row 1, col 2
#   area(2, 1, 2, 2)   # GG3 centered across both columns in row 2
# )
top_row<-(GG1 + GG2)
Second_row<-plot_spacer() + GG3 + plot_spacer()+ 
  plot_layout(widths = c(0.5, 1,0.5))
All_PCA_coords<-top_row/Second_row
All_PCA_coords
ggsave(filename = "graphiques_Surcote/theta_coords__all.png",
       plot = All_PCA_coords,
       width = 15,
       height = 10)


# Model Copula ------------------------------------------------------------

Struct_Matrice<-function_Structure_Matrice(NB_dim =ncol(resultat_obs))
Unifs<-VineCopula::pobs(resultat_obs)
chosen_Model<-VineCopula::RVineCopSelect(data =Unifs,
                                   Matrix = Struct_Matrice)

u1<-Unifs[,1]
u2<-Unifs[,2]
FAMILY_set<-c(1:40,104,114,
              124,134,204,214,
              224,234)
# Create a results data frame
FIT_cop<-function(u1,u2,fam){
  fit <- tryCatch({
    VineCopula::BiCopEst(u1, u2, family = fam)
  }, error = function(e) return(NULL))
  results<-list()
  if (!is.null(fit)) {
    k <- length(fit$par) + length(fit$par2)  # total number of parameters
    loglik <- fit$logLik
    aic <- 2 * k - 2 * loglik
    BIC<-fit$BIC
    results$FamilyName <- VineCopula::BiCopName(fam)
    results$LogLik <- loglik
    results$NumPar <- k
    results$AIC <- aic
    results$BIC<-BIC
    
  }
  return(results)
}
results_all<-sapply(FAMILY_set,FUN = FIT_cop
         ,u1=u1,u2=u2)
Val<-sapply(results_all,FUN=function(x){
  return(length(x))
})
Ind_pos<-which(Val>0)
List_pos<-lapply(Ind_pos,FUN=function(x){
  results_all[[x]]})
DF<-t(cbind(sapply(List_pos,unlist)))
DF<-as.data.frame(DF)
DF$AIC<-as.numeric(DF$AIC)
Sort_df<-DF[order(DF$AIC), ]
View(Sort_df)
write.csv(DF,file = "residus_clust/Sumy_copules.csv")

# Second 
u1<-Unifs[,1]
u2<-Unifs[,2]
FAMILY_set<-c(1:40,104,114,
              124,134,204,214,
              224,234)

results_all2<-sapply(FAMILY_set,FUN =FIT_cop,u1=Unifs[,2],
                     u2=Unifs[,3])
Val2<-sapply(results_all2,FUN=function(x){
  return(length(x))
})
Ind_pos2<-which(Val2>0)
List_pos2<-lapply(Ind_pos2,FUN=function(x){
  results_all2[[x]]})
DF2<-t(cbind(sapply(List_pos2,unlist)))
DF2<-as.data.frame(DF2)
DF2$AIC<-as.numeric(DF2$AIC)
Sort_df2<-DF2[order(DF$AIC), ]
write.csv(Sort_df2,file = "residus_clust/Sumy_copules_second.csv")

Sumy<-summary(chosen_Model)
J<-ncol(resultat_obs)
N_for_test<-nrow(Unifs)
Return_iso_density_theor_vs_obs<-function(choice_couple,J,Sumy,
                                          family_Model,Unifs,Prop_inertia,
                                          N_for_test,
                                          Orig_data){
  i<-choice_couple[1]
  j<-choice_couple[2]
  Model<-VineCopula::BiCopEst(u1 = Unifs[,i],
                       u2 = Unifs[,j],
                       family = Number)
  Simul_for_test<-VineCopula::BiCopSim(N = N_for_test,
                                       family = family_Model,
                                       par = Model$par,
                                       par2 = Model$par2)
  expr_i<-function_expression_prop_variance_j(prop_variance = Prop_inertia,
                                              j = i)
  expr_j<-function_expression_prop_variance_j(prop_variance = Prop_inertia,
                                              j = j)
  Result_MV<-Simul_for_test
  for(l in c(1,2)){
    Result_MV[,l]<-sapply(Simul_for_test[,l],function(x){
      return(quantile(Orig_data[,choice_couple[l]],x))
    })
  }
  return(Result_MV)
  # MIN_y<-min(apply(Simul_for_test,MARGIN = 1,FUN = min),
  #            apply(Unifs,MARGIN = 1,FUN = min))
  # MAX_y<-max(apply(Simul_for_test,MARGIN = 1,FUN = max),
  #            apply(Unifs,MARGIN = 1,FUN = max))
  # 
  # colnames(Simul_for_test)<-c("first","second")
  # data<-cbind.data.frame(Unifs[,i],Unifs[,j])
  # colnames(data)<-c("first","second")
  # DATA_ALL<-rbind.data.frame(data,Simul_for_test)
  # DATA_ALL$origin<-c(rep("data",nrow(data)),rep("simulations",nrow(Simul_for_test)))
  # 
  # ggplot(data=DATA_ALL,aes(x = first, y =second)) +
  #   geom_point(aes(col=origin,shape=origin))+
  #   xlim(c(-0.15,1.15))+
  #   ylim(c(-0.15,1.15))+
  #   geom_density_2d_filled(alpha = 0.5) + facet_wrap(vars(origin))+
  #   geom_density_2d(aes())+
  #   scale_color_manual(values=cols_)+
  #   scale_shape_manual(values = c("simulations"=17,"data"=19))+
  #   # scale_linetype_manual(values=c("simulations"=2,
  #   #                                "data"=1))+
  # 
  #   #ggtitle("Iso-density curves of the simulated and recorded coordinates")+
  #   #labs(linetype="Legend",fill="Level")+
  #   xlab(expr_i)+
  #   ylab(expr_j)+
  #   theme(axis.title=element_text(size=15),
  #         legend.text=element_text(size=10))+
  #   labs(col="Legend",shape="Legend",level="Level")
}
str(df_combo)

# First<-which(df_combo$Legend=="data")
# df_unif1<-as.data.frame(VineCopula::pobs(df_combo[First,c(1:3)]))
# colnames(df_unif1)<-colnames(df_combo)[c(1:3)]
# Second<-which(df_combo$Legend=="simulations")
# df_unif2<-as.data.frame(VineCopula::pobs(df_combo[Second,c(1:3)]))
# colnames(df_unif2)<-colnames(df_combo)[c(1:3)]
# DF_new<-rbind.data.frame(df_unif1,df_unif2)
# DF_new$Legend<-df_combo$Legend
ggplot(data=df_combo,aes(x=C_1,y=C_2,shape=Legend))+
  geom_point(aes(size=Legend,col=Legend))+
scale_color_manual(values=cols_)+
  geom_density_2d_filled(alpha = 0.5)+
  facet_wrap(vars(Legend))+
scale_shape_manual(values = c("simulations"=17,"data"=19))+
scale_size_manual(values=c("simulations"=0.75,"data"=1.5))+
labs(fill="Legend")
#+
#scale_size_manual(values=c("simulations"=0.75,"data"=1.5))
expr_i<-function_expression_prop_variance_j(prop_variance = Prop_inertia,
                                            j = 1)
expr_j<-function_expression_prop_variance_j(prop_variance = Prop_inertia,
                                            j = 2)
#Number<-as.numeric(VineCopula::BiCopName("Tawn270"))
Number<-2
RMV<-Return_iso_density_theor_vs_obs(choice_couple = c(1,2),
                              J = J,
                              Sumy = Sumy,
                              family_Model = Number,
                              Unifs = Unifs,Prop_inertia = Prop_inertia ,
                              N_for_test = 2000,
                              Orig_data=resultat_obs)
Min_value<-min(apply(X = resultat_obs,MARGIN = 2,FUN = min))
Max_value<-max(apply(X = resultat_obs,MARGIN = 2,FUN = max))

RankTransfo<-function(Xdata,x){
  n<-length(Xdata)
  Num<-sum(as.numeric(Xdata<x))
  Denom<-n+1
  return(Num/(Denom))
}
Param<-VineCopula::BiCopEst(u1 = Unifs[,1],u2 = Unifs[,2],
                            family = Number)
Student_copula<-copula::tCopula(param = Param$par,df = Param$par2)
xunif <- seq(0.01, 0.99, length = 100)
yunif <- seq(0.01, 0.99, length = 100)
gridcop <- expand.grid(x = xunif, y = yunif)
Dens_iso<-apply(gridcop,MARGIN = 1,
                FUN = function(x){copula::dCopula(u = x,Student_copula)})
Z_cop<-matrix(Dens_iso,nrow = length(xunif), ncol = length(yunif))
contour(xunif, yunif, Z_cop,
        nlevels = 50,
        xlab=expr_i,
        ylab=expr_j,
        col = "black")
points(Unifs[,1],Unifs[,2],
       xlab=expr_i,
       ylab=expr_j,
       pch=20,
       col="blue")

x <- seq(Min_value, Max_value, length = 100)
y <- seq(Min_value, Max_value, length = 100)
grid <- expand.grid(x = x, y = y)

KDE_dens_estimator<-function(Xi,x_){
  Objet_kde<-ks::kde(x = Xi,eval.points = x_)
  return(Objet_kde$estimate)
}
Value_found1<-sapply(grid$x,FUN = KDE_dens_estimator,Xi=resultat_obs[,1])
Value_found2<-sapply(grid$y,FUN = KDE_dens_estimator,Xi=resultat_obs[,2])
Dens_found<-Dens_iso*Value_found1*Value_found2
z <- matrix(Dens_found,nrow = length(x), ncol = length(y))
# Draw contour plot
contour(x, y, z,
      nlevels = 50,
      xlab=expr_i,
      ylab=expr_j,
      col = "black")
points(resultat_obs[,1],resultat_obs[,2],
       xlab=expr_i,
       ylab=expr_j,
      pch=20,
       col="blue")

Param2<-VineCopula::BiCopEst(u1 = Unifs[,2],
                            u2 = Unifs[,3],
                            family = Number)
Student_copula2<-copula::tCopula(param = Param2$par,
                                 df = Param2$par2)
# Draw contour plot
expr_k<-function_expression_prop_variance_j(prop_variance = Prop_inertia,
                                            j = 3)
Dens_iso2<-apply(gridcop,MARGIN = 1,
                FUN = function(x){copula::dCopula(u = x,Student_copula2)})
Z_cop2<-matrix(Dens_iso2,nrow = length(xunif), ncol = length(yunif))
contour(xunif, yunif, Z_cop2,
        nlevels = 50,
        xlab=expr_j,
        ylab=expr_k,
        col = "black")
points(Unifs[,2],Unifs[,3],
       xlab=expr_j,
       ylab=expr_k,
       pch=20,
       col="blue")

Value_found3<-sapply(grid$y,FUN = KDE_dens_estimator,Xi=resultat_obs[,3])
Dens_found2<-Dens_iso2*Value_found3*Value_found2
z2 <- matrix(Dens_found2,nrow = length(x), ncol = length(y))

contour(x, y, z2,
        nlevels = 20,
        xlab=expr_j,
        ylab=expr_k,
        col = "black")
points(resultat_obs[,2],resultat_obs[,3],
       xlab=expr_j,
       ylab=expr_k,
       pch=20,
       col="blue")

Simul_for_test<-as.data.frame(Simul_for_test)
data<-cbind.data.frame(Unifs[,1],Unifs[,2])
colnames(data)<-c("first","second")
colnames(Simul_for_test)<-colnames(data)
DATA_ALL<-rbind.data.frame(data,Simul_for_test)
DATA_ALL$origin<-c(rep("data",nrow(data)),rep("simulations",nrow(Simul_for_test)))

ggplot(data=DATA_ALL,aes(x = first, y =second)) +
geom_point(aes(col=origin,shape=origin))+
xlim(c(-0.15,1.15))+
ylim(c(-0.15,1.15))+
geom_density_2d_filled(alpha = 0.5) + facet_wrap(vars(origin))+
geom_density_2d(aes())+
scale_color_manual(values=cols_)+
scale_shape_manual(values = c("simulations"=17,"data"=19))+
# scale_linetype_manual(values=c("simulations"=2,
#                                "data"=1))+

#ggtitle("Iso-density curves of the simulated and recorded coordinates")+
#labs(linetype="Legend",fill="Level")+
theme(axis.title=element_text(size=15),
      legend.text=element_text(size=10))+
labs(col="Legend",shape="Legend",level="Level")


Model_used<-chosen_Model$family[3,1]
Model_used


# # With logit --------------------------------------------------------------
# # Tawn copula 1 and 2 -----------------------------------------------------
# 
# # = take a constant first asymp param -------------------------------------
# # + inverse of shape param of logistic. Below 1 ---------------------------
# Z<-2
# Covariates<-Sumy$edge[Z]
# First_index<-as.numeric(substr(Covariates,1,1))
# First_index
# Second_index<-as.numeric(substr(Covariates,3,3))
# Second_index
# require(VC2copula)
# ## rotations or not ?
# PAR<-Sumy$par[Z]  
# PAR2<-Sumy$par2[Z]
# c(PAR,PAR2)
# if(Sumy$family[Z]==2){
#   PAR2<-round(PAR2)
# }
# Object_copula<-VC2copula::BiCop2copula(family = Sumy$family[Z],
#                                             par =PAR ,
#                                             par2 = PAR2)
# 
# # # Asymp test with Tawn ----------------------------------------------------
# 
# VineCopula::BiCopVuongClarke(Unifs[,First_index],
#                          Unifs[,Second_index])
# Estimation_An<-"CFG"
# delta<-500
# w<-seq.int(0.01,1-(1/delta),by=1/delta)
# Theor<-copula::A(Object_copula,w)
# Emp<-copula::An.biv(Unifs[,c(First_index,Second_index)],w,
#             estimator = Estimation_An)
# Distance_CVM_data<-sum(nrow(Unifs)*(Theor-Emp)**2*w)
# M<-1000
# Bootstrap_A_An<-function(Copula_theor,n,Family_bicop,w,method_estimation){
#   Sims_<-VineCopula::pobs(copula::rCopula(n = n,copula = Copula_theor))
#   # Estimated law  ----------------------------------------------------------
#   # Copy what would have happened under the null hypothesis ----------------------------------------------------
#   Estimated_Copula_theor<-VineCopula::BiCopEst(Sims_[,1],Sims_[,2],
#                                                family=Family_bicop)
#   PAR2<-Estimated_Copula_theor$par2
#   if(Family_bicop==2){
#     PAR2<-round(PAR2)
#   }
#   in_copula_estimated<-VC2copula::BiCop2copula(Estimated_Copula_theor$family, 
#                                                Estimated_Copula_theor$par,
#                                                PAR2)
#   Theor_B<-A(in_copula_estimated,w)
#   Emp_B<-An.biv(Sims_,w,estimator = method_estimation)
#   Distance_CVM<-sum(n*(Theor_B-Emp_B)**2*w)
#   return(Distance_CVM)
# }
# Result_CVM<-replicate(M,Bootstrap_A_An(Copula_theor = Object_copula,
#                Family_bicop = Sumy$family[Z],
#                n = nrow(Unifs),w = w,method_estimation=Estimation_An))
# p_val_found<-1-mean(as.numeric(Result_CVM<Distance_CVM_data))
# plot(density(Result_CVM))
# abline(v=Distance_CVM_data,col="red")
# 
# # First two dimensions ----------------------------------------------------
# Copula_theor<-copula::pCopula(copula = Object_copula,
#                 u = Unifs[,c(First_index,Second_index)])
# 
# Copula_emp_data_eval<-copula::C.n(u=Unifs[,c(First_index,Second_index)],
#                                   X = Unifs[,c(First_index,Second_index)])
# Sn<-sum((Copula_theor-Copula_emp_data_eval)^2)
# Sn
# 
# Bootstrap_dce_emp_law<-function(Copula_theor,n,Family_bicop){
#   Sims_<-VineCopula::pobs(copula::rCopula(n = n,copula = Copula_theor))
#   # Estimated law  ----------------------------------------------------------
#   # Copy what would have happened under the null hypothesis ----------------------------------------------------
#   Estimated_Copula_theor<-VineCopula::BiCopEst(Sims_[,1],Sims_[,2],
#                                                family=Family_bicop)
#   PAR2<-Estimated_Copula_theor$par2
#   if(Family_bicop==2){
#     PAR2<-round(PAR2)
#   }
#   in_copula_estimated<-VC2copula::BiCop2copula(Estimated_Copula_theor$family, 
#                                              Estimated_Copula_theor$par,
#                                              PAR2)
#   Copula_theor_star<-copula::pCopula(copula =in_copula_estimated,
#                                 u = Sims_)
#   Cnstar<-copula::C.n(u=Sims_,X = Sims_)
#   S_star<-sum((Copula_theor_star-Cnstar)^2)
#   
#   return(S_star)
# }
# 
# 
# Set_S_star<-replicate(M,Bootstrap_dce_emp_law(Copula_theor = Object_copula,
#                                   n=nrow(Unifs),Family_bicop=Sumy$family[Z]))
# plot(density(Set_S_star))
# p_val_found<-1-mean(as.numeric(Set_S_star<Sn))
# p_val_found
# abline(v=Sn,col="red")

# copula::gofCopula(tawnT1Copula(),
#                     x =resultat_obs[,c(First_index,Second_index)],
#                     method="SnC",
#                     ties.method = "max",simulation="pb")
# gofEVCopula
# set.seed(7) # for reproducibility
# gumbel.cop <- gumbelCopula(3, dim=2)
# View(gumbel.cop)
# x <- rCopula(n, gumbel.cop) # "true" observations (simulated)
# x
# 
# u.fit <- pobs(resultat_obs[,c(First_index,Second_index)])
# Object_copula_Tawn@family
# Object<-fitCopula(tawnT1Copula(), u.fit,method="ml",
#                   estimate.variance=FALSE)
# 
# # pseudo-observations
# ## Inverting Kendall's tau
# gofT2stat(u.fit,u.fit, useR = FALSE)
# 
# object<-fitCopula( Object_copula_Tawn, Unifs)
# fcop<-object@copula
# 
# m<-500
# g <- seq(0, 1 - 1/m, by = 1/m)
# object<-copula::cramer_vonMises_Afun
# s <- .C(cramer_vonMises_Afun, as.integer(n), as.integer(m), 
#         as.double(-log(Unifs[, 1])), as.double(-log(Unifs[, 2])), as.double(A(fcop, 
#                                                                       g)), stat = double(2), as.integer(estimator == "CFG"))$stat
# return_value<-copula::A(fcop, g)
# # VineCopula::BiCopGofTest(u1 = Unifs[,First_index],
# #                          u2 = Unifs[,Second_index],
# #                          family = 134,
# #                          par = Sumy$par[1],
# #                          par2=Sumy$par2[2])
# # gofstat for testing null hypothesis  ------------------------------------
# # that it follows Unif ----------------------------------------------------
# 
# 
# without logit, if we have the students coordinates ----------------------
chosen_Model$pair.AIC
Rasymp<-VineCopula::RVineGofTest(VineCopula::pobs(resultat_obs),chosen_Model,
                       B=0)
Rasymp
# Rwhite_boostrap<-VineCopula::RVineGofTest(VineCopula::pobs(resultat_obs),chosen_Model,
#                          method="White",
#                          statistic="CvM",
#                           B = 200)
# Rwhite_boostrap

# # Berg --------------------------------------------------------------------
# # require(ADGofTest)
# RBERG<-VineCopula::RVineGofTest(VineCopula::pobs(resultat_obs),chosen_Model,
#                          method="Breymann",
#                          statistic="AD",
#                          B = 200)
# RBERG
# RECP<-VineCopula::RVineGofTest(VineCopula::pobs(resultat_obs),chosen_Model,
#                                method="ECP",
#                                B = 200)
# RECP
# # # RECP

# RECP2<-VineCopula::RVineGofTest(VineCopula::pobs(resultat_obs),chosen_Model,
#                          method="ECP2",
#                          statistic="KS",
#                          B = 200)
# # RECP2
# # RECP_2_CvM<-VineCopula::RVineGofTest(VineCopula::pobs(resultat_obs),chosen_Model,
# #                          method="ECP2",
# #                          statistic="CvM",
# #                          B = 200)
# 



# BR vs acp ---------------------------------------------------------------
H_Params<-c("radial","logit",50)
NB_times<-100
ALPHA_PROP<-0.10
K<-4
typeS<-"simple"
resultat<-matrix(NA,ncol=3,nrow=3)

source("fonctions/fonctions_perfs_ML.R")
# Use of the function -----------------------------------------------------
Lperfs_ACP<-Running_perfs_ML(Base_simul = simulations,Base_data = obs_exts, 
                         hyp_param = H_Params, NB_times =100,alpha_prop=ALPHA_PROP,K = K, 
                         type_sampling = typeS,title_ROC =paste0(nom_present," ACP"),
                         NB_shown_ROC = 20)

# BR ----------------------------------------------------------------------
Lperfs_BR<-Running_perfs_ML(Base_simul =simulations_BR,Base_data = obs_exts, 
                         hyp_param = H_Params, NB_times =100,alpha_prop=ALPHA_PROP,K = K, 
                         type_sampling = typeS,title_ROC =paste0(nom_present," BR"),
                         NB_shown_ROC = 20)
Lperfs_ACP$Accuracy
Lperfs_BR$Accuracy
Norm<-c(apply(X = simulations,MARGIN = 1,FUN = calcul_norme_L2),
                        apply(X = simulations,MARGIN = 1,FUN = calcul_norme_L2),
                        apply(X =obs_exts,MARGIN = 1,FUN = calcul_norme_L2))
Type<-c(rep("ACP",nrow(simulations)),
      rep("BR",nrow(simulations_BR))                             ,
      rep("data",nrow(obs_exts)))
Type
df_norm_type<-cbind.data.frame(Norm,Type)
colnames(df_norm_type)<-c("norm","type")
ggplot(data=df_norm_type,aes(y=norm,x=type,fill=type))+
geom_violin(alpha=0.75)+
geom_boxplot(width=0.25)+
scale_fill_manual(values=c("data"="blue",
                           "ACP"="orange",
                           "BR"="cyan"))+
labs(fill="Legend")+
ylab("L2 norm (m)")

# Same density for BR/ACP -------------------------------------------------
# so it does not affect the ML differences ----------------------------

# ANGLE -------------------------------------------------------------------
Angle_ACP<-t(t(simulations)%*%diag(apply(simulations,FUN = calcul_norme_L2,
                                   MARGIN=1)^(-1)))
Angle_BR<-t(t(simulations_BR)%*%diag(apply(simulations_BR,FUN = calcul_norme_L2,
                                   MARGIN=1)^(-1)))
Angle_exts<-t(t(obs_exts)%*%diag(apply(obs_exts,FUN = calcul_norme_L2,
                                     MARGIN=1)^(-1)))
Lperfs_A_ACP<-Running_perfs_ML(Base_simul = Angle_ACP,Base_data = Angle_exts, 
                           hyp_param = H_Params, NB_times =100,alpha_prop=ALPHA_PROP,K = K, 
                           type_sampling = typeS,title_ROC =paste0(nom_present," ACP"),
                           NB_shown_ROC = 20)
Lperfs_A_BR<-Running_perfs_ML(Base_simul =Angle_BR,Base_data = Angle_exts, 
                          hyp_param = H_Params, NB_times =100,alpha_prop=ALPHA_PROP,K = K, 
                          type_sampling = typeS,title_ROC =paste0(nom_present," BR"),
                          NB_shown_ROC = 20)
Lperfs_A_ACP$Accuracy
Lperfs_A_BR$Accuracy

# ACP rules ! ---------------------------------------------------------------

# # Fin du code -------------------------------------------------------------
# #----------------------------------------------------------------------
# 
# data<-cbind.data.frame(Unifs[,1],Unifs[,2])
# colnames(data)<-c("first","second")
# Unifs_really_used<-cbind.data.frame(VineCopula::pobs(resultat_simul[,c(1:2)]))
# colnames(Unifs_really_used)<-colnames(data)
# DATA_ALL<-rbind.data.frame(data,
#                            Unifs_really_used)
# DATA_ALL$origin<-c(rep("data",nrow(data)),rep("simulations",nrow(Unifs_really_used)))
# 
# ggplot(data=DATA_ALL,aes(x = first, y =second)) +
#   xlim(c(-0.15,1.15))+
#   ylim(c(-0.15,1.15))+
#   geom_density_2d_filled(alpha = 0.5) + facet_wrap(vars(origin))+
#   geom_density_2d(aes(linetype=origin))+
#   scale_linetype_manual(values=c("simulations"=2,
#                                  "data"=1))+
#   
#   #ggtitle("Iso-density curves of the simulated and recorded coordinates")+
#   labs(linetype="Legend",fill="Level")+
#   theme(axis.title=element_text(size=15))

