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
require(tibble)
require(threshr)

### Source des fonctions. ####
source("fonctions.R")
F_names<-ls()
nb_coeurs<-detectCores()-4
coeurs<-makeCluster(nb_coeurs)
clusterExport(coeurs,varlist=F_names)

fonction_stationnarité<-function(coeurs,nom_variable,opt_tend,liste_t,type_donnees,lims_y_gamma,debut_P,pas){
  repertoire<-"ss_tend/"
  if(type_donnees=="HIVER"){
    repertoire<-paste0(repertoire,"HIVER/")
  }
  if(type_donnees=="HORS_HIVER"){
    repertoire<-paste0(repertoire,"HORS_HIVER/")
  }
  nom_recherche<-paste0(repertoire,nom_variable,"_ss_tend.csv")
  obs<-read.csv(nom_recherche)[,2:38]
  Tendance_ounon<-ifelse(opt_tend,"On enlève la tendance","On conserve la tendance")
  print(Tendance_ounon)
  if(opt_tend==FALSE){
    #Nous avons séparé la tendance lineaire et l'aspect aleatoire. 
    #Si nous négligeons la tendance, nous travaillons sur les 
    #données brutes. 
    #L'addition est donc necessaire. 
    nom_rech2<-paste0(repertoire,"RL_",nom_variable,".csv")
    RL_result<-read.csv(nom_rech2)[,2:38]
    obs<-obs+RL_result
  }
  colnames(obs)<-c(1:length(colnames(obs)))
  print(nrow(obs))
  boxplot(obs)
  VECTEUR_normeL2<-parApply(cl=coeurs,X = obs,MARGIN = 1,FUN = calcul_norme_L2)
  ACF_norme<-acf(VECTEUR_normeL2,plot=FALSE,lag=1000)
  PACF_norme<-pacf(VECTEUR_normeL2,plot=FALSE,lag=1000)
  plot(ACF_norme,ylim=c(-1,1),main=paste("ACF de la norme de",nom_variable,"pour",type_donnees))
  plot(PACF_norme,ylim=c(-1,1),main=paste("PACF de la norme de",nom_variable,"post tri",type_donnees))
  for (t in liste_t){
    variable_t<-unlist(obs[,t])
    df_time_variable<-cbind.data.frame(1:length(variable_t),variable_t)
    colnames(df_time_variable)<-c("time","obs")
    ACF_t<-acf(variable_t,lag=1000,plot = FALSE)
    PACF_t<-pacf(variable_t,lag=1000,plot = FALSE)
    vecty<-(30:length(variable_t)/2)
    plot(ACF_t,ylim=c(-1,1),main=paste("ACF de",nom_variable,"en t=",t,"post tri, pour",type_donnees))
    plot(PACF_t,ylim=c(-1,1),main=paste("PACF de",nom_variable,"en t=",t,"post tri, pour",type_donnees))
    #results_ml<-parSapply(cl=coeurs,vecty,fonction_ML_extRemes,variable_t)
    #print(fonction_MLplot_resume(resultatML = results_ml,vecteur_k = vecty,nom_variable = paste0(nom_variable,"en t= ",t),lims_Y=lims_y_gamma))
  }
  # Analyse des tendances des données ---------------------------------------
  resultat<-parallel::parApply(cl=coeurs,X = obs,MARGIN = 2,FUN = FNCT_percentile_ppourcent,pas=pas,debut=debut_P)
  Noms_lignes<-c("Moyenne","Ecart_type")
  lim<-nrow(resultat)-3
  for(l in c(1:lim)){
    Noms_lignes<-c(Noms_lignes,paste0(as.character((pas*l+debut_P)*100),"%"))
  }
  Noms_lignes<-c(Noms_lignes,"99%")
  df<-t(as.data.frame(resultat))
  colnames(df)<-Noms_lignes
  M<-melt(df)
  colnames(M)<-c("Temps","variable","valeur")
  M1<-M[which((M$variable!="Ecart_type")),]
  print(M1)
  GG<-ggplot(data = M1,aes(x=Temps,y=valeur,group=interaction(variable),col=variable))+
    geom_point()+
    geom_line()+
    ggtitle(paste0("Description des données pour ",nom_variable," en hiver"))
  print(GG)
  M2<-M[which(M$variable%in%c("Moyenne")),]
  print(M2)
  GG1<-ggplot(data = M2,aes(x=Temps,y=valeur))+
    geom_point()+
    geom_line()+
    ggtitle(paste0("Moyenne des données pour ",nom_variable," en hiver"))
  print(GG1)
  return(NA)
}
a<-fonction_stationnarité(coeurs=coeurs,nom_variable = "Hs",opt_tend = TRUE,liste_t = 1,type_donnees = "HIVER",lims_y_gamma=c(-2,2),debut_P = 0.90,pas=0.025)
#fonction_stationnarité(coeurs=coeurs,nom_variable = "U",opt_tend = TRUE,liste_t = 1,type_donnees = "HIVER")
#fonction_stationnarité(coeurs=coeurs,nom_variable = "Surcote",opt_tend = TRUE,liste_t = 1,type_donnees = "HIVER")

stopCluster(coeurs)