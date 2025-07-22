rm(list=ls())
set.seed(133)

#### Packages. #####
require(corrplot)
require(ggplot2)
require(dplyr)
require(reshape2)
require(FactoMineR)
require(fda)
require(extRemes)
require(parallel)
require(POT)

### Source des fonctions. ####
source("fonctions/fonctions_r_pareto_EVOLTAU.R")
F_names<-ls()
nb_coeurs<-detectCores()-4
coeurs<-makeCluster(nb_coeurs)
clusterExport(coeurs,varlist=F_names)
annees<-seq.int(1979,2016,by=1)
dataFiles<-c()

for(elt in annees){
  commentaire<-paste0("donnees_transformees/evt_",elt,"*.csv")
  file_annee <- parLapply(cl=coeurs,Sys.glob(commentaire), read.csv)
  dataFiles<-c(dataFiles,file_annee)
}


# Extraction Johanna ------------------------------------------------------
dates_J<-c("08/03/2008","09/03/2008","10/03/2008")
creation_df_sans_tendance<-function(df,nom_variable,coeurs,SELECT,saut_pas,dates_Johanna){
  obs<-parLapply(cl=coeurs,df,f_nom_variable,nom_variable)
  
  # ACF de la serie des normes (donnees brutes) ----------------------------------------------
  Svariable<-parSapply(cl=coeurs,obs,calcul_norme_L2)
  ACF_S<-acf(Svariable,lag=1000,plot = FALSE)
  PACF_S<-pacf(Svariable,lag=1000,plot = FALSE)
  plot(ACF_S,ylim=c(-1,1),main=bquote("ACF of N"[.(nom_variable)]))
  plot(PACF_S,ylim=c(-1,1),main=bquote("PACF of N"[.(nom_variable)]))
  
  # ACF en chaque temps.  ---------------------------------------------------
  for (t in 1:1){
    fnct_t<-function(coord,obs){
      return(obs[coord])
    }
    variable_t<-sapply(obs,fnct_t,coord=t)
    ACF_t<-acf(variable_t,lag=400,plot = FALSE)
    PACF_t<-pacf(variable_t,lag=400,plot = FALSE)

    # TSA::periodogram(variable_t,plot=TRUE,xlab="Fréquence",ylab="Périodogramme")
    # abline(v=1/730,col="green")
    #title(paste0("Analyse de la variable ",nom_variable," en t=",t))
    vecty<-(30:length(variable_t)/2)
    plot(ACF_t,ylim=c(-1,1),main="",cex.lab=1.5)
    #paste("ACF of",nom_variable,"at t=",t)
    plot(PACF_t,ylim=c(-1,1),main="",cex.lab=1.5)
    #main=paste("PACF of",nom_variable,"at t=",t)
  }
  TEMPS<-parLapply(cl=coeurs,df,f_nom_variable,"Date")
  intermed<-t(do.call(cbind.data.frame,TEMPS))
  rownames(intermed)<-1:nrow(intermed)

  # Date pic de la marée.  --------------------------------------------------
  date_compar<-intermed[,19]
  colonne<-as.numeric(substr(date_compar,4,5))
  indices_J<-which(date_compar%in% dates_Johanna)
  DATE_Johanna<-date_compar[indices_J]
  # Recuperer obs Johanna ---------------------------------------------------
  
  liste_Johanna<-cbind.data.frame(obs[indices_J])
  n_rep<-ncol(liste_Johanna)
  matin_soir<-rep(c("(morning)","(evening)"),n_rep)
  date_matin_soir<-paste(DATE_Johanna,matin_soir)
  Rt<-t(liste_Johanna)
  rownames(Rt)<-c(1:6)
  resultat<-melt(Rt)
  colnames(resultat)<-c("indu","temps","valeur")
  resultat$date<-date_matin_soir[resultat$indu]
  GGJohanna_brut<-ggplot(data=resultat,aes(x=temps,y=valeur,col=date,group=interaction(date)))+
    geom_point()+
    geom_line()+
    ggtitle(paste0("Observations durant la tempête Johanna pour ",nom_variable))
  print(GGJohanna_brut)
  vecteur_hiver<-c(1:3,9:12)
  Indices_autres<-which(!colonne%in%vecteur_hiver)
  Indices_hiver<-which(colonne%in%vecteur_hiver)
  date_compar_hiver<-date_compar[Indices_hiver]
  indices_J<-which(Indices_hiver%in%indices_J)
  obs_hiver<-obs[Indices_hiver]

  # export winter dates -----------------------------------------------------
  write.csv(file = "ss_tend/HIVER/dates.csv",x=date_compar_hiver)
  
  #ACF obtenu
  variable_hiver1<-sapply(obs_hiver,fnct_t,coord=1)
  ACF_t<-acf(variable_hiver1,lag=1000,plot = FALSE)
  plot(ACF_t,ylim=c(-1,1),main=paste("ACF of",nom_variable,"at t=1 (winter)"))
  
  L_hiver<-length(obs_hiver)
  obs_autres<-obs[Indices_autres]
  # Hiver -------------------------------------------------------------------
  # regression lineaire -----------------------------------------------------
  vecteur_reg_lin<-c()
  vecteur_r<-c()
  vecteur_signf<-c()
  vecteur_beta_chap<-c()
  vecteur_r_J<-c()
  for (col in c(1:37)){
    
    # Hiver -------------------------------------------------------------------
    variable<-sapply(obs_hiver,FUN = function(t,serie){return(serie[t])},t=col)
    L<-length(variable)
    vecteur_temps_milieu<-c(1:L)
    reg<-lm(variable~vecteur_temps_milieu)
    objet<-coef(summary(reg))
    Constante<-as.numeric(predict(object = reg,newdata=data.frame(vecteur_temps_milieu=0)))
    predictions_RL_hiver<-as.numeric(reg$fitted.values)
    Residus<-variable-predictions_RL_hiver+Constante
    vecteur_r<-c(vecteur_r,Residus)
    Residus_J<-Residus[indices_J]
    vecteur_r_J<-c(vecteur_r_J,Residus_J)
    vecteur_reg_lin<-c(vecteur_reg_lin,predictions_RL_hiver)
    p_valeur_signf<-objet[2,4]
    vecteur_beta_chap<-c(vecteur_beta_chap,objet[2,1])
    vecteur_signf<-c(vecteur_signf,p_valeur_signf)
  }
  df_coeff<-cbind.data.frame(vecteur_beta_chap)
  colnames(df_coeff)<-c(paste("Coefficient lineaire de",nom_variable,"(Hiver)"))
  df_signf<-cbind.data.frame(vecteur_signf)
  colnames(df_signf)<-c(paste("p_valeur de",nom_variable,"(Hiver)"))
  Residus_hiver<-matrix(vecteur_r,nrow=L_hiver)
  Residus_Johanna<-matrix(vecteur_r_J,nrow=length(indices_J))

  df_date_J<-cbind.data.frame(DATE_Johanna,Residus_Johanna)
  matrice_tendance<-matrix(vecteur_reg_lin,nrow=L_hiver)
  write.csv(df_date_J,file=paste0("ss_tend/Johanna_detrend_",nom_variable,".csv"))
  # Johanna -----------------------------------------------------------------
  
  Tempete_Johanna<-read.csv(paste0("ss_tend/Johanna_detrend_",nom_variable,".csv"))
  fin<-ncol(Residus_hiver)+2
  colnames(Tempete_Johanna)[3:39]<-1:ncol(Residus_hiver)
  Tempete_Johanna<-Tempete_Johanna[,2:fin]
  colnames(Tempete_Johanna)[1]<-"date"
  date<-Tempete_Johanna$date
  n_rep<-length(date)/2
  matin_soir<-rep(c("(matin)","(soir)"),n_rep)
  date_matin_soir<-paste(date,matin_soir)
  MELT<-melt(t(Tempete_Johanna[,2:38])) 
  colnames(MELT)<-c("temps","indu_stat","valeur")
  MELT$indu_stat<-date_matin_soir[MELT$indu_stat]
  GGJ<-ggplot(data=MELT,aes(x=temps,y=valeur,col=indu_stat,group=interaction(indu_stat)))+
    geom_point()+
    geom_line()+
    labs(color="Date")+
    ggtitle(paste0("Observations lors de la tempête Johanna sans la tendance pour ", nom_variable))
  print(GGJ)
  TEST_moyenne<-t(cbind.data.frame(obs_hiver))
  for (t in 1){

    # Général -----------------------------------------------------------------
    variable_t<-sapply(obs,fnct_t,coord=t)
    ACF_t<-acf(variable_t,lag=1000,plot = FALSE)
    PACF_t<-pacf(variable_t,lag=1000,plot = FALSE)
    vecty<-(30:length(variable_t)/2)
    if(nom_variable=="Surcote"){
      plot(ACF_t,ylim=c(-1,1),main="")
      plot(PACF_t,ylim=c(-1,1),main="")

    }
    else{
      plot(ACF_t,ylim=c(-1,1),main=paste("ACF of",nom_variable,"at t=",t," (whole data)"))
      plot(PACF_t,ylim=c(-1,1),main=paste("PACF of",nom_variable,"at t=",t," (whole data)"))
      
    }
    # Hiver -----------------------------------------------------------------
    variable_t<-Residus_hiver[,t]
    ACF_t<-acf(variable_t,lag=1000,plot = FALSE)
    PACF_t<-pacf(variable_t,lag=1000,plot = FALSE)
    vecty<-(30:length(variable_t)/2)
    plot(ACF_t,ylim=c(-1,1),main=paste("ACF of",nom_variable," detrended at t=",t," (hiver)"))
    plot(PACF_t,ylim=c(-1,1),main=paste("PACF of",nom_variable,"detrended at t=",t," (hiver)"))
    
  }
  
  # Saut à effectuer. -------------------------------------------------------
  Indice_pris_hiver<-seq.int(from = 1,to = L_hiver,by = saut_pas)
  vect_matin_soir<-rep(c("(morning)","(evening)"),length(date_compar_hiver)/2)
  date_compar_hiver_MSOIR<-paste(date_compar_hiver,vect_matin_soir)
  date_compar_final<-date_compar_hiver_MSOIR[Indice_pris_hiver]
  Residus_detrend_hiver<-Residus_hiver[Indice_pris_hiver,]
  Mat<-matrice_tendance[Indice_pris_hiver,]
  NL2<-apply(Residus_detrend_hiver,MARGIN = 1,FUN = calcul_norme_L2)
  indice<-which.max(NL2)
  write.csv(x=date_compar_final,file="ss_tend/HIVER/dates_prises.csv")
  
  # ACF post selection --------------------------------------------------------
  Nvariable<-parApply(cl=coeurs,Residus_detrend_hiver,calcul_norme_L2,MARGIN=1)
  ACF_S<-acf(Nvariable,lag=1000,plot = FALSE)
  PACF_S<-pacf(Nvariable,lag=1000,plot = FALSE)
  plot(ACF_S,ylim=c(-1,1),main=bquote("ACF après transformation de N"[.(nom_variable)]))
  ACF_residus<-acf(Residus_detrend_hiver[,1],lag=1000,plot = FALSE)
  PACF_residus<-pacf(Residus_detrend_hiver[,1],lag=1000,plot = FALSE)
  plot(ACF_residus,ylim=c(-1,1),main=paste("ACF of",nom_variable," detrended at t=",t," (winter), one out of",saut_pas))
  print(length(Nvariable))
  plot(PACF_residus,ylim=c(-1,1),main=paste("PACF of",nom_variable," detrended at t=",t," (winter), one out of",saut_pas))
  
  # Presentation résultat ---------------------------------------------------

  # (1) Stat desc ---------------------------------------------------------------
  mu_t<-colMeans(Residus_detrend_hiver)
  sigma_t<-apply(X = Residus_detrend_hiver,MARGIN = 2,FUN = sd)
  df_stat_desc<-cbind.data.frame(mu_t,sigma_t)
  colnames(df_stat_desc)<-c("Moyenne","Ecart-type")
  df_stat_desc<-melt(t(df_stat_desc))
  colnames(df_stat_desc)<-c("Statistique","Temps","Valeur")
  TITRESD<-paste0("Statistiques pour ",nom_variable, " en hiver")
  STATdesc<-ggplot(data = df_stat_desc,aes(x=Temps,y=Valeur,col=Statistique,group=interaction(Statistique)))+
    geom_line()+
    geom_point()+
    ggtitle(TITRESD)
    
  print(STATdesc)
  
  B<-melt(t(Residus_detrend_hiver))
  colnames(B)<-c("Temps","Individu","Mesure")
  B$Temps<-as.factor(B$Temps)
  TITRE<-paste0("Densité en chaque temps pour ",nom_variable, " en hiver")
  PV<-ggplot2::ggplot(data = B,aes(Temps,Mesure,fill=Temps))+
    geom_violin()+
    geom_boxplot(width=0.2)+
    guides(fill = "none")+
    ggtitle(TITRE)+
    ylab("Densité")
  print(PV)  
  plot(colMeans(TEST_moyenne),main=paste0("Fonction moyenne pour ",nom_variable," en hiver"))
  plot(colMeans(Residus_detrend_hiver),main=paste0("Fonction moyenne pour les résidus de ",nom_variable," en hiver"))
  write.csv(x = Residus_detrend_hiver,file = paste0("ss_tend/HIVER/",nom_variable,"_ss_tend.csv"))
  write.csv(x=Mat,file=paste0("ss_tend/HIVER/RL_",nom_variable,".csv"))
  
  return(list("Signf"=df_signf,"Coeff"=df_coeff,"Johanna"=Rt,"dates_matin_soir"=date_matin_soir))
}
#saut de 6 normalement. 
# pas de 3 pour alterner marée du matin et marée de l'après-midi. 
# essai avec pas de 1.
# R<-creation_df_sans_tendance(df=dataFiles,nom_variable = "U",coeurs=coeurs,SELECT=TRUE,
#                              saut_pas =3,dates_Johanna = dates_J)
# 
# R2<-creation_df_sans_tendance(df=dataFiles,nom_variable = "Hs",coeurs=coeurs,SELECT=TRUE,
#                               saut_pas=3,dates_Johanna = dates_J)

R3<-creation_df_sans_tendance(df=dataFiles,nom_variable = "Surcote",coeurs=coeurs,SELECT=TRUE,
                              saut_pas = 3,dates_Johanna = dates_J)


melteur_<-rbind(R$Johanna,R2$Johanna,R3$Johanna)
Nb_rep<-nrow(melteur_)/3

test<-melt(t(melteur_))
test$variable<-c(rep("U",37*6),rep("Hs",37*6),rep("Surcote",37*6))

colnames(test)<-c("Temps","individu","valeur","variable")
test$individu<-R$dates_matin_soir[test$individu]

ggplot(data=test,aes(x=Temps,y = valeur,linetype=variable,col=individu,group=interaction(individu,variable)))+
  geom_point()+
  geom_line()+
  facet_wrap(~variable,scales="free_y")+
  ggtitle(paste0("Observations obtenues lors de la tempête Johanna sans la tendance"))
### Fermeture parallel ------------------------------------------------------
stopCluster(coeurs)

#
Tendance_ttes_vars<-cbind.data.frame(R$Signf,R2$Signf,R3$Signf)
resultat<-melt(t(Tendance_ttes_vars))
colnames(resultat)<-c("Variable","Temps","p_valeur")
GGsignf<-ggplot(data = resultat,aes(x=Temps,y=p_valeur,col=Variable,group=interaction(Variable)))+
  geom_point()+
  geom_line()+
  ylab("p valeur")+
  xlab("Temps")+
  geom_hline(yintercept = 0.05,show.legend = TRUE)+
  ggtitle("p valeur du test (Ho): aucune tendance linéaire (1979-2016)")+
  labs(caption = "p valeur du test à 5%")

print(GGsignf)
Tendance_ttes_beta<-cbind.data.frame(R$Coeff,R2$Coeff,R3$Coeff)
resultat<-melt(t(Tendance_ttes_beta))
colnames(resultat)<-c("Variable","Temps","beta")

GGbeta<-ggplot(data = resultat,aes(x=Temps,y=beta,col=Variable,group=interaction(Variable)))+
  geom_point()+
  geom_line()+
  ylab("beta")+
  xlab("Temps")+
  ggtitle("Coefficient linéaire estime (1979-2016)")
print(GGbeta)

