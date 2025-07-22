rm(list=ls())
set.seed(133)
par(mfrow=c(2,2))
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
  
fonction_par_tempsSSTEND<-function(coeurs,type_seuils,nom_variable,opt_tend,liste_t,lMin_u,lMax_u,NCLUST,U,type_donnees){
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
  boxplot(obs)
  VECTEUR_normeL2<-parApply(cl=coeurs,X = obs,MARGIN = 1,FUN = calcul_norme_L2)
  vecteur_resultat<-list()
  vecteur_pu<-list()
  j<-1
  for (t in liste_t){
    #On veut pour chaque temps donne en entree approcher la quantite
    #mu permettant d'approcher la loi des excedents par une Pareto.
    variable_t<-unlist(obs[,t])
    df_time_variable<-cbind.data.frame(1:length(variable_t),variable_t)
    colnames(df_time_variable)<-c("time","obs")
    if(NCLUST!=0){
      #Les correlations entre les observations peuvent
      #compliquer la tache.
      #La methode clust de POT permet de prendre en
      #compte un decalage.
      print(clust(data=df_time_variable,u =U,tim.cond = NCLUST,plot = FALSE,clust.max =FALSE))
      variable_t<-clust(data=df_time_variable,u =U,tim.cond = NCLUST,plot = FALSE,clust.max = TRUE)[,"obs"]
      df_time_variable<-cbind.data.frame(1:length(variable_t),variable_t)
      colnames(df_time_variable)<-c("time","obs")
    }
    # #Plusieurs criteres sont utilises:
    # #(a) tcplot. Estimation des parametres d'une GPareto dont le seuil est mu1.
    # # Il faut que les parametres estimes soient constants 
    # #pour tout mu1>mu. 
    # #(b) mrlplot. Graphique de l'esperance residuelle en fonction de mu1. 
    # #La fonction de mu1 doit etre lineaire lorsque mu1>mu. 
    # #(c) diplot. Le nombre d'excedents au cours du temps
    # #peut etre approche par un processus de Poisson si mu est bien choisi.
    # #On compare le parametre de la loi avec la moyenne empirique. 
    # #Le rapport doit etre proche de 1.
    if(type_seuils=="Threshr"){
      Min_u<-quantile(variable_t,lMin_u[t])
      Max_u<-quantile(variable_t,lMax_u[t])
      qu_threshr<- quantile(variable_t, probs = seq(lMin_u[t], lMax_u[t],length.out=10))
    }
    if(type_seuils=="POT"){
      Min_u<-lMin_u[t]
      Max_u<-lMax_u[t]
      qu_threshr<-seq.int(from=Min_u,to=Max_u,length.out=10)
    }
    
    par(mfrow=c(2,2))
    valeurs<-POT::tcplot(variable_t,ask = FALSE,u.range = c(Min_u,Max_u))
    titre<-paste0("t=",t," ",nom_variable," sans la tendance en ",tolower(type_donnees))
    title(titre,line = -3, outer = TRUE)
    POT::mrlplot(data = variable_t,u.range = c(Min_u,Max_u))
    # # With time.

    POT::diplot(data = df_time_variable,u.range = c(Min_u,Max_u),nt=100)
    Resultat_threshr<-threshr::ithresh(data=variable_t,u_vec =qu_threshr)
    par(mfrow=c(1,1))
    vecteur_pu[[j]]<-summary(Resultat_threshr)[4]
    vecteur_resultat[[j]]<-summary(Resultat_threshr)[3]
    j<-j+1
  }
  plot(unlist(vecteur_pu))
  vecteur_resultat<-unlist(vecteur_resultat)
  vect_SEUILS<-(100-unlist(vecteur_pu))/100
  
  ss_actend<-ifelse(opt_tend," sans tendance"," avec tendance")
  plot(vect_SEUILS,main=paste0("poids (p_u) mis en chaque temps pour ",nom_variable,ss_actend),type="o",ylab="Poids",xlab="Temps")
  Mean_pu<-mean(vect_SEUILS)
  Seuils_finaux<-rep(NA,length(vect_SEUILS))
  #return(vecteur_resultat)
  if(type_seuils=="POT"){
    return(vecteur_resultat)
  }
  else{
    for (t in liste_t){
      variable_t<-unlist(obs[,t])
      Seuil_t<-quantile(variable_t,1-Mean_pu)
      Seuils_finaux[t]<-Seuil_t
    }
    return(Seuils_finaux)
  }
  
}
fnct_lancement<-function(vecteur_min,vecteur_max,nom_variable,NCLUST,U,type_seuils){
  SEUILS_THRESHR<-fonction_par_tempsSSTEND(nom_variable=nom_variable,opt_tend=TRUE,coeurs=coeurs,liste_t=c(1:37),lMin_u=vecteur_min,lMax_u=vecteur_max,NCLUST=NCLUST,U=U,type_donnees = "HIVER",type_seuils=type_seuils)
  write.table(SEUILS_THRESHR,file = paste0("SEUILS_THRESHR_",nom_variable,"_HIVER_",type_seuils,".txt"))
}
# HIVER -------------------------------------------------------------------
#vecteur_min_h<-rep(0.85,37)
#vecteur_max_h<-rep(0.98,37)
# 4 et 5.5 sinon
# 4.6
vecteur_min_h<-rep(4.6,37)
vecteur_max_h<-rep(4.9,37)
#Il faut au moins 50 individus.
fnct_lancement(vecteur_max = vecteur_max_h,vecteur_min = vecteur_min_h,nom_variable = "Hs",NCLUST=0,U=2,
               type_seuils="POT")

# U ---------------------------------------------------------------
vecteur_q1<-rep(0.85,37)
vecteur_q2<-rep(0.98,37)
# # vecteur_minh<-rep(12,37)
# # vecteur_maxh<-rep(14,37)
# # # Essayer avec (12,14) pour les seuils. 
#fnct_lancement(vecteur_max = vecteur_q2,vecteur_min = vecteur_q1,nom_variable = "U",NCLUST=0,U=0,type_seuils="Threshr")

# Surcote -----------------------------------------------------------------
# 0.2,0.26
vecteur_Q1<-rep(0.85,37)
vecteur_Q2<-rep(0.96,37)
# # vecteur_Q1<-rep(0.2,37)
# # vecteur_Q2<-rep(0.26,37)
#fnct_lancement(vecteur_max = vecteur_Q2,vecteur_min = vecteur_Q1,nom_variable = "Surcote",NCLUST = 0,U = 0,type_seuils="Threshr")
# # # #


# Suppression des coeurs --------------------------------------------------
stopCluster(coeurs)
