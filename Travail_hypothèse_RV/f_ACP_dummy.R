Simul_dummy<-function(coeurs,p_norm,nom,opt_tend,BORNES_Y=c(0,30),type_donnees,rho_hissage){
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
  colnames(obs)<-c(1:length(colnames(obs)))
  vecteur_norm_L2<-apply(obs,FUN = calcul_norme_L2,MARGIN = 1)
  Q<-quantile(vecteur_norm_L2,1-p_norm)
  Indices_select<-which(vecteur_norm_L2>Q)
  OBS_EXT<-obs[Indices_select,]
  VECT_EXT<-vecteur_norm_L2[Indices_select]
  individus<-t(do.call(cbind.data.frame,OBS_EXT))
  FORME_d<-t(individus%*%diag(VECT_EXT^(-1)))
  
  fnct_centrer_par_colonne<-function(x){
    return(x-mean(x))
  }
  mu_t<-colMeans(FORME_d)
  observations_centrees<-apply(X = FORME_d,FUN = fnct_centrer_par_colonne,MARGIN = 2)
  nlignes<-nrow(observations_centrees)
  Mat_moyenne<-matrix(rep(mu_t,nlignes),byrow = TRUE,nrow=nlignes)
  sigma_ANGLE<-apply(FORME_d,MARGIN = 2,FUN = sd)
  mu_ANGLE<-apply(FORME_d,MARGIN = 2,FUN = mean)
  
  ACP_inds_extremes<-PCA(X = observations_centrees,graph = FALSE,scale.unit =FALSE)
  F_propres_normal<-ACP_inds_extremes$svd$V
  Coordonnees<-ACP_inds_extremes$ind$coord
  val_lambda<-ACP_inds_extremes$eig[,1]
  vect_diff<-cumsum(c(0,val_lambda))/sum(val_lambda)
  vecteur_propvarexp<-1-vect_diff
  plot(vecteur_propvarexp,main=paste0("Evolution de la contribution des axes scale=TRUE ",nom," en ",type_donnees),ylab="Proportion non expliquee",type="o")
  
  Scores_1_base<-Coordonnees[,1]
  Scores_2_base<-Coordonnees[,2]
  Coordonnees_ACP_normal<-cbind(Scores_1_base,Scores_2_base)
  fpropre_base<-F_propres_normal[,1]
  fpropre_base2<-F_propres_normal[,2]
  plot(density(Scores_1_base))
  plot(density(Scores_2_base))
  # 
  M<-nrow(OBS_EXT)
  valeurs_unif<-runif(M)
  reals_P_std<-sapply(valeurs_unif,Generation_Pareto_std)
  
  NB_OBS<-length(Scores_1_base)
  modele<-fitdistrplus::fitdist(Scores_1_base,"norm")
  plot(modele)
  # modele_st<-fitdistrplus::fitdist(Scores_1_base,"t",start=)
  # plot(modele_st)
  Scores_simules_Premdim<-rnorm(NB_OBS,mean=mean(Scores_1_base),sd = sd(Scores_1_base))
  plot(density(Scores_simules_Premdim))
  Matrice_score<-diag(Scores_simules_Premdim)
  Matrice_repmfonction<-matrix(rep(fpropre_base,M),ncol=M)
  Simulations<-t(Matrice_repmfonction%*%Matrice_score)
  FORME_ACP<-Simulations+Mat_moyenne
  
  # Predictions_simul -------------------------------------------------------
  Pred_ACP<-t(t(FORME_ACP)%*%diag(reals_P_std))
  par(mfrow=c(1,2))
  matplot(individus,type="l")
  NOM<-paste0("Individus sélectionnées pour ",nom," en ",tolower(type_donnees))
  title(NOM)
  matplot(t(Pred_ACP),type="l")
  NOM2<-paste0("Simulations (ACP) pour ",nom)
  title(NOM2)
  par(mfrow=c(1,1))

  # # Alternative -------------------------------------------------------------
  # ACP_inds_extremes<-PCA(X = t(individus),graph = FALSE,scale.unit =TRUE)
  # F_propres_normal<-ACP_inds_extremes$svd$V
  # Coordonnees<-ACP_inds_extremes$ind$coord
  # val_lambda<-ACP_inds_extremes$eig[,1]
  # vect_diff<-cumsum(c(0,val_lambda))/sum(val_lambda)
  # vecteur_propvarexp<-1-vect_diff
  # plot(vecteur_propvarexp,main=paste0("Evolution de la contribution des axes scale=FALSE ",nom," en ",type_donnees),ylab="Proportion non expliquee",type="o")
  # 
  # Scores_1_base<-Coordonnees[,1]
  # Scores_2_base<-Coordonnees[,2]
  # Coordonnees_ACP_normal<-cbind(Scores_1_base,Scores_2_base)
  # fpropre_base<-F_propres_normal[,1]
  # fpropre_base2<-F_propres_normal[,2]
  # plot(density(Scores_1_base))
  # plot(density(Scores_2_base))
}
