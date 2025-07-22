fonction_export_convergence<-function(coeurs,nom,n.dens,opt_tend,Tau,BORNES_Y=c(0,30),type_donnees){
  
  # Import des donnees ------------------------------------------------------
  repertoire<-"ss_tend/"
  if(type_donnees=="HIVER"){
    repertoire<-paste0(repertoire,"HIVER/")
  }
  if(type_donnees=="HORS_HIVER"){
    repertoire<-paste0(repertoire,"HORS_HIVER/")
  }
  # Tau_en _t ---------------------------------------------------------------)
  nom_rech1<-paste0(repertoire,nom,"_ss_tend.csv")
  obs<-read.csv(nom_rech1)
  colnames(obs)<-c(1:38)
  obs<-obs[,2:38]
  ## Table pour la tendance --------------------------------------------------
  nom_rech2<-paste0(repertoire,"RL_",nom,".csv")
  MATRICE_tend<-read.csv(nom_rech2)[,2:38]
  if (opt_tend==FALSE){
    obs<-obs+MATRICE_tend
  }
  else{
    observations_AVECTEND<-obs+MATRICE_tend
    colnames(observations_AVECTEND)<-c(1:length(colnames(observations_AVECTEND)))
  }
  p_u<-c()
  ## Approcher le poids mis dans la queue de distribution.
  for(j in c(1:ncol(obs))){
    prop<-mean(as.numeric(obs[,j]>Tau[j]))
    p_u<-c(p_u,prop)
  }
  colnames(obs)<-c(1:length(colnames(obs)))
  ## Preciser si la tendance est utilisee.
  Vecteur_DIM<-1:ncol(obs)
  
  # Comparaison marginales Pareto, modele de melange ------------------------
  Resultat_P<-as.data.frame(t(parSapply(cl=coeurs,Vecteur_DIM,f_marginales_all_Pareto,donnees=obs,p_u=p_u,n.dens=n.dens,donnees_a_transformer=obs)))
  K<-Resultat_P$kernel_dens
  PARETO<-Resultat_P$obs
  vecteur_pareto<-do.call(cbind.data.frame,PARETO)
  rownames(vecteur_pareto)<-1:nrow(vecteur_pareto)
  colnames(vecteur_pareto)<-c(1:ncol(vecteur_pareto))
  
  # Preparation des donnees pour ACP  --------------------------------------------------------------------
  ###################################################
  fnct_centrer_par_colonne<-function(x){
    return(x-mean(x))
  }
  OBS<-vecteur_pareto
  ### Travail RV pour les scores.
  Resultat<-FactoMineR::PCA(OBS,scale.unit = FALSE,ncp = 2,graph = FALSE)
  Scores_acp<-Resultat$ind$coord
  ###
  Coord<-Scores_acp[,1]
  vect_y<-20:250
  ML_extRemes<-parSapply(cl = coeurs,X =vect_y ,FUN = fonction_ML_extRemes,donnees=Coord)
  LIMS<-c(-1,3)
  print(fonction_MLplot_resume(resultatML = ML_extRemes,vecteur_k = vect_y,nom_variable = paste0("le score C1 de T(",nom,")"),lims_Y=LIMS))
  
  Coord2<-Scores_acp[,2]
  ML_extRemes<-parSapply(cl = coeurs,X =vect_y ,FUN = fonction_ML_extRemes,donnees=Coord2)
  print(fonction_MLplot_resume(resultatML = ML_extRemes,vecteur_k = vect_y,nom_variable = paste0("le score C2 de T(",nom,")"),lims_Y=LIMS))
  write.csv(OBS,file="OBS_PARETO.csv")
}

