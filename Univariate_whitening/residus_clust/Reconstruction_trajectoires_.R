rm(list=ls())
set.seed(133)
par(mfrow=c(1,1))
source("../fonctions/fonctions.R")
require(reshape2)
require(ggplot2)
require(gridExtra)
require(grid)
require(dplyr)
require(extRemes)
require(ggside)
require(scales)

cols_<-c("data"="blue","simulations"="orange","confidence_band"="darkblue")
methode_<-"fenetre"
#Obtenir à partir des simulations  ---------------------------------------
# des résidus de nouvelles trajectoires -----------------------------------
Nom_v<-"Surcote"
inds_extremes<-paste0("inds_extremes_donnees_",Nom_v,".csv")

Reconstruction_trajectoire<-function(fichier_Simulation_residus,nom_variable,fich_reg1,type_donnees,fichier_idx_extremes, 
                                     dates_Johanna){
  nom_present<-ifelse(Nom_v=="Surcote","Surge",Nom_v)
  repertoire<-"../ss_tend/"
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
  obs<-read.csv(nom_recherche)[,2:38]
  colnames(obs)<-c(1:37)
  
  # Import des individus extrêmes dans la base ------------------------------
  fichier_idx_extremes<-paste0("residus_clust/",fichier_idx_extremes)
  Inds_extremes_L2<-read.csv(fichier_idx_extremes)[,2]
  Individus_extremes_base<-obs[Inds_extremes_L2,]
  rownames(Individus_extremes_base)<-c(1:nrow(Individus_extremes_base))
  NL2_residus_exts<-NL2_residus[Inds_extremes_L2]

  # Import dates ------------------------------------------------------------
  Dates_prises<-read.csv(file="../ss_tend/HIVER/dates_prises.csv")[,2]
  Dates_prises_extremes<-Dates_prises[Inds_extremes_L2]
  mini_date<-substr(Dates_prises_extremes,1,10)
  indice_Joh<-which(mini_date%in%dates_Johanna)
  date_J_found<-Dates_prises_extremes[indice_Joh]
  # On prend la dernière obs ------------------------------------------------
  Simulations_epsi<-read.csv(file=fichier_Simulation_residus)[,2:38]
  NL2_epsi_simul<-apply(X=Simulations_epsi,
                        MARGIN=1,FUN = calcul_norme_L2)
  par(mfrow=c(1,2))
  colnames(Simulations_epsi)<-c(1:37)
  matplot(t(Simulations_epsi),type="l")
  matplot(t(Residus[Inds_extremes_L2,]),type="l")
  par(mfrow=c(1,1))

  # Eventual relation between  epsilon_t et residus -------------------------------------------
  number_obs<-nrow(obs)-1
  df_resultat<-matrix(NA,nrow = 37,ncol = 4)
  DF_EV_resid<-read.csv(file=paste0("../Travail_hypothèse_RV/EVA_",Nom_v,"_residus.csv"))[,2:5]
  vect_threshold_GPD<-DF_EV_resid$seuil_t
  
  # Correlation_t=19 --------------------------------------------------------
  val_peak<-obs[1:number_obs,19]
  r_peak<-Residus[c(2:nrow(obs)),19]
  
  for(t in c(1:ncol(obs))){

    # subset of extreme/Mann Whitney ------------------------------------------
    seuil_t<-vect_threshold_GPD[t]
    series_r<-Residus[c(2:nrow(obs)),t]
    veille_r<-obs[1:number_obs,t]
    indices_ext<-which(series_r>seuil_t)
    indice_non_ext<-which(series_r<=seuil_t)
    subset1<-veille_r[indices_ext]
    subset2<-veille_r[indice_non_ext]

    # Correlation test --------------------------------------------------------

    Pearson_t<-cor.test(x = series_r,y=veille_r,method="pearson")$p.value
    Spearman_t<-cor.test(x = series_r,y=veille_r,method="spearman")$p.value
    Kendall_t<-cor.test(x = series_r,y=veille_r,method="kendall")$p.value
    MaWhitney_U<-wilcox.test(subset1,subset2)$p.value
    df_resultat[t,]<-c(Pearson_t, Spearman_t,
                       Kendall_t,MaWhitney_U)
  }
  df_resultat<-as.data.frame(df_resultat)
  colnames(df_resultat)<-c("Pearson_t","Spearman_t", 
                           "Kendall_t","Mann_Whitney_t")
  
  df_melt_corr<-melt(t(df_resultat))
  colnames(df_melt_corr)<-c("method","Time","value")
  titre<-expression("Correlations between "~espilon[M]~" and "~X[M-1]~"for"~Nom_v)
  GG_corr<-ggplot(data=df_melt_corr,aes(x=Time,y=value,group=interaction(method),col=method))+
    geom_line()+
    geom_point()+
    xlab("Time (hour) with respect to tidal peak")+
    ylab("p value")+
    geom_hline(yintercept = 0.05,col="black")+
    ggtitle(label = do.call("substitute", list(titre[[1]], list(Nom_v = Nom_v))))
   
  ## replace the variable by its value with the do.call and the list. 
  print(GG_corr)
  
  # Modele AR ---------------------------------------------------------------
  Mod_AR<-read.csv(file=fich_reg1)
  Fin<-ncol(Mod_AR)
  Mod_AR<-Mod_AR[,c(2:Fin)]
  
  
  M<-nrow(Simulations_epsi)
  Y_Mat<-matrix(NA,nrow = M,ncol = 37)
  echantillonnage_YO<-function(niv_simul,echs_N2){
    return(which.min(abs(niv_simul-echs_N2)))
  }
  # Se baser sur des obs non exts -------------------------------------------
  Indices_non_ext<-which(!c(1:nrow(obs))%in%Inds_extremes_L2==TRUE)
  
  NL2<-apply(X = obs,MARGIN = 1,FUN = calcul_norme_L2)
  inds_non_nuls<-which(NL2_residus>0)
  plot(NL2_residus[inds_non_nuls],NL2[inds_non_nuls],
       main="Relation des deux normes", 
       xlab=expression("Norme de epsilon"[M]),ylab=expression("Norme de Y"[M]))
  y<-NL2[inds_non_nuls]
  x<-NL2_residus[inds_non_nuls]
  modele_lin_epsi_Y<-lm(y~x)
  JBtest<-tseries::jarque.bera.test(modele_lin_epsi_Y$residuals)
  seuil_y<-quantile(y,0.75)
  seuil_x<-quantile(x,0.75)

  # fonction_craft ---------------------------------------------------------
  ##################
  number_obs<-length(NL2)-1
  NL2_veille<-NL2[1:number_obs]
  # plot(density(NL2),main="densité norme L2")
  filtre<-min(NL2[Inds_extremes_L2])

  prob_non_exceedance<-mean(as.numeric(NL2<=filtre))
  print("Return period")
  npy<-length(NL2)/37
  periode_retour<-(npy*(1-prob_non_exceedance))^(-1)
  Indices_non_ext<-which(NL2<filtre)
  par(mfrow=c(1,1))
  N_XO<-c()

  # Choix des indices -------------------------------------------------------
  # Ici on prend les observations ayant une norme inférieure au minimum des normes
  # des observations extrêmes (norme résidu). 
  df_relation_dsbase<-cbind.data.frame(NL2_residus[2:length(NL2_residus)],
                                       NL2_veille)
  
  colnames(df_relation_dsbase)<-c("epsi","L2_indu_pris")
  if(methode_=="non_ext")
  {
    Echantillons_pris<-sample(x =Indices_non_ext,
                              size = M,replace=TRUE)
  }
  if(methode_=="free"){
    Echantillons_pris<-sample(x =c(1:nrow(obs)),
                              size = M,replace=TRUE)
  }
  inds_chosen<-which((df_relation_dsbase$epsi>0))
  df_subset<-df_relation_dsbase[inds_chosen,]
  if(methode_=="emp_copula"){
    df_emp_C<-cbind(NL2_veille,NL2_residus[2:length(NL2_residus)])
    couple_unif<-VineCopula::pobs(df_emp_C)
    EmpCopula<-copula::empCopula(couple_unif)
    return(df_emp_C)
  }
  if(methode_=="copule"){
    couple_unif<-VineCopula::pobs(df_subset)
    Family_couple_L2<-VineCopula::BiCopSelect(couple_unif[,1],couple_unif[,2])
    print("########### Famille choisie #############")
    print(summary(Family_couple_L2))
    
    qsimul<-sapply(NL2_epsi_simul,FUN=function(x){
      return(mean(as.numeric(x>df_subset$epsi)))
    })
    fnct_simul<-function(x_simul,nl2_Yo){
      u<-VineCopula::BiCopCondSim(1, cond.val = x_simul, cond.var = 1, Family_couple_L2)
      return(quantile(nl2_Yo,u))
    }
    value_predicted<-sapply(X = qsimul,fnct_simul,
                            nl2_Yo=df_subset$L2_indu_pris)
    Echantillons_pris<-sapply(X = value_predicted,echantillonnage_YO,
                              echs_N2=NL2)
  }
  if(methode_=="fenetre"){
    ech_fenetre<-sapply(X =NL2_epsi_simul,Sample_window,
           y=NL2_veille,x_window=NL2_residus[2:length(NL2_residus)],size_window=20)
    Echantillons_pris<-ech_fenetre
    
  }
  if(methode_=="H_Tawn"){
    # Travail Gaussian estimation ---------------------------------------------
    seuil_epsi<-quantile(df_subset$epsi,0.75)
    seuil_Y<-quantile(df_subset$L2_indu_pris,0.75)
    value_predicted<-fonction_simul_HTawn(x=df_subset$epsi,y=df_subset$L2_indu_pris, 
                         seuil_x = seuil_epsi,seuil_y = seuil_Y, 
                         vecteur_x_reg =NL2_epsi_simul)
    Echantillons_pris<-sapply(X = value_predicted,echantillonnage_YO,
                              echs_N2=NL2)
  }
  
  # alternative : trouver avec modele lm -------------------------------------
  L<-nrow(obs)-1
  sub_obs<-NL2[1:L]
  sub_resid<-NL2_residus[2:nrow(obs)]
  sub_obs<-sub_obs[which(sub_resid>0)]
  sub_resid<-subset(sub_resid,sub_resid>0)
  df<-cbind.data.frame(sub_obs,sub_resid)
  colnames(df)<-c("obs_passe","resid")
  modele_choix<-lm(obs_passe~resid,data = df)
  predictions_niveaux_passe<-predict(modele_choix,data.frame(resid=NL2_epsi_simul))

  par(mfrow=c(1,2))
  plot(density(NL2[Echantillons_pris]),main="norme L2 des obs (X0) prises")
  plot(density(NL2[Inds_extremes_L2-1]),main="norme L2 veille des obs extremes")
  par(mfrow=c(1,1))
  nb_decal<-0
  if(ncol(Mod_AR)>2){
    nb_decal<-ncol(Mod_AR)-2
  }
  Simulate_extreme_from_resid<-function(nb_decal,Sample_taken,
                                        obs,Mod_AR,j,Simul_epsi){
    indice<-sort(c(Sample_taken[j])+c(-nb_decal:0),decreasing = TRUE)
    individu_pris<-obs[indice,]
    Prediction<-as.numeric(sapply(X = c(1:ncol(obs)),FUN = Function_AR_p,
                                  obs_eve=individu_pris,
                                  epsilon_t_plus1=unlist(Simul_epsi[j,]),
                                  Modele_ar_par_t = Mod_AR))
    return(Prediction)
  }
  Y_Mat<-t(sapply(c(1:M),Simulate_extreme_from_resid,nb_decal = nb_decal,
                       Sample_taken = Echantillons_pris ,
                       Mod_AR = Mod_AR,Simul_epsi = Simulations_epsi,
                       obs = obs))
  Y_Mat<-as.data.frame(Y_Mat)
  colnames(Y_Mat)<-c(1:ncol(Y_Mat))
  NL2_simulations<-apply(X =  Y_Mat,MARGIN = 1,FUN=calcul_norme_L2)
  print("simul vs veille--YO")
  test<-Inds_extremes_L2[2]-1
  # Simulation evenement type Johanna ---------------------------------------
  Veille_Johann<-obs[Inds_extremes_L2[indice_Joh]-1,]
  date_veille_J<-Dates_prises[Inds_extremes_L2[indice_Joh]-1]
  Johann<-obs[Inds_extremes_L2[indice_Joh],]
  M_Prediction_obs_J<-t(sapply(c(1:nrow(Simulations_epsi)),
                              FUN = function(x)
      {vecteur<-as.numeric(sapply(X = c(1:ncol(obs)),FUN = Function_AR_p,
                                                     obs_eve=Veille_Johann,
                                                     epsilon_t_plus1=unlist(Simulations_epsi[x,]),
                                                     Modele_ar_par_t = Mod_AR))
       return(vecteur)}))
  
  matplot(t(M_Prediction_obs_J),type="l")
  mtext(paste0("Prediction pour Johanna"))
  lines(x=c(1:37),Johann,col="red",lwd=3,lty=2)
  
  lendemain<-obs[Inds_extremes_L2[indice_Joh]+1,]

  M_Prediction_jour_extreme<-t(sapply(c(1:nrow(Simulations_epsi)),
                             FUN = function(x)
                             {vecteur<-as.numeric(sapply(X = c(1:ncol(obs)),FUN = Function_AR_p,
                                                         obs_eve=Johann,
                                                         epsilon_t_plus1=unlist(Simulations_epsi[x,]),
                                                         Modele_ar_par_t = Mod_AR))
                             return(vecteur)}))
  matplot(t(M_Prediction_jour_extreme),type="l")
  mtext(paste0("Prediction pour le ",date_J_found," pour ",Nom_v))
  
  # versus observation normal -----------------------------------------------
  
  X_zero_nomal<-obs[Echantillons_pris[1],]
  date_XOnormal<-Dates_prises[Echantillons_pris[1]]
  lendemain_normal<-obs[Echantillons_pris[1]+1,]
  M_Prediction_jour_normal<-t(sapply(c(1:nrow(Simulations_epsi)),
           FUN = function(x){vecteur<-as.numeric(sapply(X = c(1:ncol(obs)),FUN = Function_AR_p,
                                       obs_eve=X_zero_nomal,
                                       epsilon_t_plus1=unlist(Simulations_epsi[x,]),
                                       Modele_ar_par_t = Mod_AR))
                 return(vecteur)}))
  
  matplot(t(M_Prediction_jour_normal),type="l")
  mtext("preds_jour_normal")
  
  # Relation veille_XO vs norme epsilon -------------------------------------
  NL2_epsi_simul<-apply(X = Simulations_epsi,MARGIN = 1,
                        FUN = calcul_norme_L2)
  df_residus_vs_veille_simul<-cbind.data.frame(NL2_epsi_simul,NL2[Echantillons_pris])
  colnames(df_residus_vs_veille_simul)<-c("simul_epsi","L2_indu_pris")
  df_residus_vs_veille_simul$indicatrice<-sapply(Echantillons_pris,function(x){return(x%in%Inds_extremes_L2)})

    # Simulations -------------------------------------------------------------
  compar_pic<-cbind.data.frame(r_peak,val_peak)
  colnames(compar_pic)<-c("r_peak","val_peak")
  GG_rel_compar<-ggplot(data = df_relation_dsbase,aes(y=L2_indu_pris,
                                                                 x=epsi))+
    geom_point()+
    geom_point(data=compar_pic,aes(x=r_peak,y=val_peak),pch=2)+
    ggtitle("Compar link norm vs peak BW epsi and previous day")

  print(GG_rel_compar)
  GG_rel_v2<-ggplot(data=df_residus_vs_veille_simul,aes(x=simul_epsi,
                                                        y=L2_indu_pris,col="simulations"))+
    geom_point(pch=17)+
    xlab(expression("L2 norm of "~epsilon[M]))+
    ylab(expression("L2 norm of "~tilde(X)[M-1]))+
    geom_point(data = df_relation_dsbase,aes(y=L2_indu_pris,
                                             x=epsi,col="data"))+
    scale_color_manual(values=c("simulations"="orange","data"="blue"))+
    labs(col="Legend")+
    #ggtitle(expression("Representation of the couple ("~epsilon[M]~","~tilde(X)[M-1]~")"))
    theme(axis.title=element_text(size=15),
          legend.text=element_text(size=10))

  print(GG_rel_v2)
  
  # Comparaison des extremogrammes ------------------------------------------
  Matrice_couples<-list()
  L<-ncol(Y_Mat)
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
  Tau1<-sapply(c(1:37),function(x,q,t){
    return(as.numeric(quantile(x[,t],q)))},q=q_chosen,
    x=Y_Mat)
  Tau2<-sapply(c(1:37),function(x,q,t){
    return(as.numeric(quantile(x[,t],q)))},q=q_chosen,
    x=Individus_extremes_base)
  
  # Extremogramme donnees est proche de 1 car mauvais choix de Tau. 
  # extrapolation donc le Tau doit être aussi extrapolé. 
  
  Extremogram_simulations<-empirical_extremogram(Matrix_couples = Matrice_couples,
                                                 inds_select = Y_Mat,
                                                 Tau = Tau1)
  Extremogram_donnees<-empirical_extremogram(Matrix_couples = Matrice_couples,inds_select= Individus_extremes_base,
                                             Tau = Tau2)
  N_rep<-500
  B<-nrow(Individus_extremes_base)
  Resultat<-replicate(n =N_rep,expr = fnct_estim_extremo_resample(B = B,inds_ext  = Individus_extremes_base,
                                                               Tau = Tau2,Matrix_couples = Matrice_couples,
                                                               vector_distances=vecteur_distances))

  estimateur_95<-apply(Individus_extremes_base,MARGIN =2,function(x){return(as.numeric(quantile(x,0.95)))})
  estimateur_975<-apply(Individus_extremes_base,MARGIN =2,function(x){return(as.numeric(quantile(x,0.975)))})
  estimateur_05<-apply(Individus_extremes_base,MARGIN =2,function(x){return(as.numeric(quantile(x,0.05)))})
  estimateur_50<-apply(Individus_extremes_base,MARGIN =2,function(x){return(as.numeric(quantile(x,0.50)))})
  
  # estimateur simulations --------------------------------------------------
  estimateur_95_simul<-apply(Y_Mat,MARGIN =2,function(x){return(as.numeric(quantile(x,0.95)))})
  estimateur_05_simul<-apply(Y_Mat,MARGIN =2,function(x){return(as.numeric(quantile(x,0.05)))})
  estimateur_50_simul<-apply(Y_Mat,MARGIN =2,function(x){return(as.numeric(quantile(x,0.50)))})
  estimateur_975_simul<-apply(Y_Mat,MARGIN =2,function(x){return(as.numeric(quantile(x,0.975)))})
  L_q<-c(0.95,0.975,0.05,0.50)
  Tendencies<-replicate(n =N_rep,expr = Resamples_tendencies(B = B,extremes_indus = Individus_extremes_base,
                                                      list_Q = L_q))
  borne_plus95<-apply(Tendencies[1,,],MARGIN = 1,FUN = function(x){return(as.numeric(quantile(x,0.975)))})
  borne_moins95<-apply(Tendencies[1,,],MARGIN = 1,FUN = function(x){return(as.numeric(quantile(x,0.025)))})

  borne_plus975<-apply(Tendencies[2,,],MARGIN = 1,FUN = function(x){return(as.numeric(quantile(x,0.975)))})
  borne_moins975<-apply(Tendencies[2,,],MARGIN = 1,FUN = function(x){return(as.numeric(quantile(x,0.025)))})
  
  borne_plus05<-apply(Tendencies[3,,],MARGIN = 1,FUN = function(x){return(as.numeric(quantile(x,0.975)))})
  borne_moins05<-apply(Tendencies[3,,],MARGIN = 1,FUN = function(x){return(as.numeric(quantile(x,0.025)))})

  borne_plus50<-apply(Tendencies[4,,],MARGIN = 1,FUN = function(x){return(as.numeric(quantile(x,0.975)))})
  borne_moins50<-apply(Tendencies[4,,],MARGIN = 1,FUN = function(x){return(as.numeric(quantile(x,0.025)))})
  df<-cbind.data.frame(c(estimateur_05,estimateur_50,estimateur_95,estimateur_975),
                       c(borne_moins05,borne_moins50,
                         borne_moins95,borne_moins975),
                       c(borne_plus05,borne_plus50,
                         borne_plus95,borne_plus975))
  colnames(df)<-c("estimateur","borne_moins","borne_plus")

  df$estimateur_simul<-c(estimateur_05_simul,estimateur_50_simul
                         ,estimateur_95_simul,estimateur_975_simul)
  df$Temps<-rep(c(1:ncol(Individus_extremes_base)),4)
  df$Temps<-(df$Temps-19)*(1/6)
  df$percent<-c(rep("Q05",ncol(Individus_extremes_base)),
                 rep("Q50",ncol(Individus_extremes_base)),
                 rep("Q95",ncol(Individus_extremes_base)), 
                 rep("Q975",ncol(Individus_extremes_base)))
  GG_percent<-ggplot(data=df,aes(x=Temps,y=estimateur,group=interaction(percent),col="data"))+
    facet_wrap(~percent,scales="free_y")+
    geom_line()+
    geom_point()+
    geom_ribbon(mapping = aes(ymin=borne_moins,ymax=borne_plus,col="confidence_band"),alpha=0.15,
                fill="grey", linetype = "dashed")+
    geom_line(aes(x=Temps,y=estimateur_simul,col="simulations"),linetype=2)+
    #ggtitle(paste0("Comparison between extremes observations and observations, ",nom_present))+
    theme(axis.title=element_text(size=15),
          legend.text=element_text(size=10))+
    ylab("Surge (m)")+
    scale_color_manual(values=cols_)+
    xlab("Time (hour) with respect to tidal peak")+
    labs(col="Legend")
  print(GG_percent)
    
    
  fnct_quantile<-function(x,q){
    M<-apply(X = x,MARGIN = 2,FUN = function(X){return(as.numeric(quantile(X,q)))})
    return(M)
  }
  Q05<-fnct_quantile(x = t(Resultat),q=0.025)
  Q95<-fnct_quantile(x = t(Resultat),q=0.975)
  df<-cbind.data.frame(Extremogram_simulations,Extremogram_donnees,vecteur_distances)
  colnames(df)<-c("valeur_simulation","valeur_donnees","delta")
  
  resultat_delta<-df %>% group_by(delta) %>% summarise(val_donnees=mean(valeur_donnees),
                                                       val_simul=mean(valeur_simulation))
  resultat_delta$borne_inf<-Q05
  resultat_delta$borne_sup<-Q95
  travail_extremo<-as.data.frame(resultat_delta[,c("val_donnees","val_simul","borne_inf","borne_sup")])
  travail_extremo$Temps<-c(1:nrow(travail_extremo))
  GG_extremo<-ggplot(travail_extremo,aes(x=Temps,y=val_donnees,col="data"))+
    geom_line()+
    geom_point()+
    geom_line(aes(y=val_simul,col="simulations"))+
    geom_point(aes(y=val_simul,col="simulations"),pch=2)+
    geom_ribbon(mapping = aes(ymin=borne_inf,ymax=borne_sup,col="confidence_band"),alpha=0.15,
                fill="grey", linetype = "dashed")+
    scale_color_manual(values=cols_)+
    #ggtitle(paste0("Estimated extremogram for ",nom_present," with resampling"))+
    theme(axis.title=element_text(size=15),
          legend.text=element_text(size=10))+
    labs(col="Legend")+
    xlab("Lag h")+
    ylab("value")
  
  print(GG_extremo)
  L_predJ<-list("simulations"=M_Prediction_obs_J,
                "X0"=Veille_Johann, 
                "date_X0"=date_veille_J, 
                "X1"=Johann, 
                "date_Johanna"=date_J_found)
  lJ<-list("simulations"=M_Prediction_jour_extreme,
           "X0"=Johann, 
           "date_J"=date_J_found, 
           "X1"=lendemain)
  lnormal<-list("simulations"=M_Prediction_jour_normal,
                "X0"=X_zero_nomal, 
                "X1"=lendemain_normal, 
                "date_XO"=date_XOnormal)
  return(list("simulations_AR"=Y_Mat,"observations_exts_base"=Individus_extremes_base,
              "dates_evenements_extremes"=Dates_prises_extremes,
              "liste_Johanna"=lJ,
              "liste_normal"=lnormal,
              "X"=obs,
              "indices_X0"=Echantillons_pris, 
              "indices_exts"=Inds_extremes_L2, 
              "liste_Pred_Johanna"=L_predJ))

}
dates_Johanna<-c("10/03/2008","08/03/2008", 
                 "09/03/2008")
methode<-"ACP"
if(methode=="ACP"){
  fichier_simul<-paste0("residus_clust/Simulations",Nom_v,"_avec_std.csv")
}
if(methode=="BR")
{  
  fichier_simul<-paste0("residus_clust/Simulations",Nom_v,"_avec_std_BR.csv")
}

nom_present<-ifelse(Nom_v=="Surcote","Surge",Nom_v)
fich_reg<-paste0("differenciation_travail/",Nom_v,"_Reg_AR1.csv")
resultat_Simul_Obs<-Reconstruction_trajectoire(fichier_Simulation_residus = fichier_simul,
                           nom_variable = Nom_v,
                           type_donnees="HIVER",
                           fich_reg1 =fich_reg,
                           fichier_idx_extremes=inds_extremes, 
                           dates_Johanna=dates_Johanna)

# Construction apprentissage-test -----------------------------------------
Simulations<-resultat_Simul_Obs$simulations_AR
matplot(t(Simulations),type="l")
Realite<-resultat_Simul_Obs$observations_exts_base
mu1<-apply(X = Simulations,MARGIN = 2,FUN = mean)
mu2<-apply(X = Realite,MARGIN = 2,FUN = mean)


Nb_ech<-nrow(Realite)
Nb_simul<-nrow(Simulations)
Base_ind_s<-c(1:Nb_simul)
Base_ind_r<-c(1:Nb_ech)
Indices_train_r<-sample(Base_ind_r,size = round(Nb_ech*0.7),replace = FALSE)
Indices_test_r<-which(!(Base_ind_r%in% Indices_train_r))
length(Indices_test_r)
Indices_s<-sample(Base_ind_s,size =Nb_ech,replace = FALSE)
nb_train<-round(0.70*Nb_ech)
Indices_train_s<-Indices_s[c(1:nb_train)]
debut<-nb_train+1
Indices_test_s<-Indices_s[c(debut:Nb_ech)]

# travail choix Xo --------------------------------------------------------
simul_Xoest<-resultat_Simul_Obs$liste_Johanna$simulations[Indices_test_s,]
Joh<-as.numeric(resultat_Simul_Obs$liste_Johanna$X0)
X1<-as.numeric(resultat_Simul_Obs$liste_Johanna$X1)
date_X0ext<-resultat_Simul_Obs$liste_Johanna$date_J

simul_Xoest<-melt(t(simul_Xoest))
ncol(simul_Xoest)
colnames(simul_Xoest)<-c("time","id","value")
simul_Xoest$id<-as.character(simul_Xoest$id)
simul_Xoest$X_zero<-Joh[simul_Xoest$time]
simul_Xoest$X_un<-X1[simul_Xoest$time]
simul_Xoest$time<-(simul_Xoest$time-19)/6
GGJoh<-ggplot(data=simul_Xoest,aes(x=time,y=value,
                            group=interaction(id), 
                            col=id))+
  geom_line(alpha=0.7)+
  guides(col="none")+
  xlab(" ")+ylab(" ")+
  geom_line(aes(y=X_zero),col="red", 
            linetype="dotted",size=1.2)+
  labs(linetype="Date of X0")+
  annotate("text", x = 0, y = X1[30], 
           label=date_X0ext)
GGJoh
norm<-as.numeric(resultat_Simul_Obs$liste_normal$X0)
X1_normal<-as.numeric(resultat_Simul_Obs$liste_normal$X1)
simul_Xoest_normal<-melt(t(resultat_Simul_Obs$liste_normal$simulations[Indices_test_s,]))
colnames(simul_Xoest_normal)<-c("time","id","value")
simul_Xoest_normal$id<-as.character(simul_Xoest_normal$id)
simul_Xoest_normal$X_zero<-norm[simul_Xoest_normal$time]
simul_Xoest_normal$X_un<-X1_normal[simul_Xoest_normal$time]
simul_Xoest_normal$time<-(simul_Xoest_normal$time-19)/6
Date_normal<-resultat_Simul_Obs$liste_normal$date_XO

GG_normal<-ggplot(data=simul_Xoest_normal,aes(x=time,y=value,
                            group=interaction(id), 
                            col=id))+
  geom_line(alpha=0.7,show.legend=FALSE)+
  xlab(" ")+ylab(" ")+
  geom_line(aes(y=X_zero),col="red", 
            linetype="dotted",size=1.2)+
  annotate("text", x = 0, y = X1_normal[30], 
           label=Date_normal)

min_grid<-min(min(simul_Xoest$value), min(simul_Xoest_normal$value), 
              min(simul_Xoest$X_zero),min(simul_Xoest_normal$X_zero), 
              X1_normal[30],X1[30])
max_grid<-max(max(simul_Xoest$value), max(simul_Xoest_normal$value), 
              max(simul_Xoest$X_zero),max(simul_Xoest_normal$X_zero),
              X1_normal[30],X1[30])

GG_normal<-GG_normal+ylim(c(min_grid,max_grid))
GGJoh<-GGJoh+ylim(c(min_grid,max_grid))
# combined <-GGJoh + GG_normal & theme(legend.position = "bottom")
# combined + plot_layout(guides = "collect")
yleft <- textGrob(paste0(nom_present," (m)"),
                  rot=90,
                  gp = gpar(col = "black", fontsize = 15))
Xbottom<-textGrob("Time (hour) with respect to tidal peak",
                  gp = gpar(col = "black", fontsize = 15))

gridExtra::grid.arrange(GGJoh,GG_normal,ncol=2,
                        left=yleft,bottom=Xbottom)

Indices_exts<-read.csv(paste0("residus_clust/",inds_extremes))[,2]
Indices_exts
# Analyse des résultats ---------------------------------------------------

Simul<-resultat_Simul_Obs$simulations_AR

NL2_simul<-apply(X = Simul,MARGIN = 1,FUN = calcul_norme_L2)
Obs_ext<-resultat_Simul_Obs$observations_exts_base
NL2_ext<-apply(X = Obs_ext,MARGIN = 1,FUN = calcul_norme_L2)
matplot(t(Simul),type="l")
inds_non_nuls<-which(NL2_simul>0)
print(length(inds_non_nuls))
Angle_inds<-t(t(Obs_ext)%*%diag(NL2_ext^(-1)))
Angle_simul<-t(t(Simul[inds_non_nuls,])%*%diag(NL2_simul[inds_non_nuls]^(-1)))

# PCA pour analyser coords ------------------------------------------------

# direct PCA --------------------------------------------------------------

PCA_brut_<-FactoMineR::PCA(X = Obs_ext,scale.unit = TRUE,graph = TRUE)

sd<-apply(X =Obs_ext,MARGIN = 2,FUN = sd)
mu<-colMeans(Obs_ext)
resultat_brut<-scale(Simul,center = mu,scale = sd)
conversion_donnees_brut<-resultat_brut%*%PCA_brut_$svd$V[,c(1:2)]
coords_donnees_brut<-PCA_brut_$ind$coord[,c(1:2)]

percent_variance_brut<-PCA_brut_$eig[c(1:2),2]

ks.test(conversion_donnees_brut[,1],coords_donnees_brut[,1])
wilcox.test(conversion_donnees_brut[,1],coords_donnees_brut[,1])
plot(coords_donnees_brut)
points(conversion_donnees_brut,col="red")

# Reject KS test if we use directly PCA -----------------------------------
# Reject null hypothesis for the first dimension (KS and MW) --------------
ks.test(conversion_donnees_brut[,2],coords_donnees_brut[,2])
wilcox.test(conversion_donnees_brut[,2],coords_donnees_brut[,2])

df_ACP_simul<-as.data.frame(conversion_donnees_brut)
colnames(df_ACP_simul)<-sapply(c(1:ncol(df_ACP_simul)), 
                               function(x){return(paste0("Score_",x))})

df_ACP_data<-as.data.frame(coords_donnees_brut[Indices_test_r,])
colnames(df_ACP_data)<-sapply(c(1:ncol(df_ACP_data)), 
                               function(x){return(paste0("Score_",x))})

# define a new df to use a categorical variable ---------------------------

df_combo<-rbind.data.frame(df_ACP_data,df_ACP_simul)
df_combo$type<-c(rep("data",nrow(df_ACP_data)),rep("simulations",nrow(df_ACP_simul)))
GG0<-ggplot(data=df_ACP_data,aes(x=Score_1,y=Score_2,col="data"),size=1.5,pch=19)+
  geom_point()+
  geom_point(data=df_ACP_simul,aes(x=Score_1,y=Score_2,col="simulations"),
             size=0.75,
             pch=17)+
  geom_xsidedensity(data=df_combo,aes(fill=type), alpha = 0.5)+
  geom_ysidedensity(data=df_combo,aes(fill=type), alpha = 0.5)+
  xlab(paste0("First dimension (",round(percent_variance_brut[1],1),"% of the variance)"))+
  ylab(paste0("Second dimension (",round(percent_variance_brut[2],1),"% of the variance)"))+
  labs(col="Legend")+
  scale_color_manual(values=cols_)+
  scale_fill_manual(values=cols_)+
  guides(fill="none")+
  theme(axis.title=element_text(size=15),
        legend.text=element_text(size=10))
  
GG0
#theme(text = element_text(size=rel(3.5)))
#ggtitle(paste0("Coordinates in the PCA basis of ",ifelse(Nom_v=="Surcote","S",Nom_v)," for a sample (n=",length(Indices_test_r),")"))

# Keep for the null hypothesis for the second dimension BUT little variance -----------------------

# PCA on the angles -------------------------------------------------------
PCA_inds<-FactoMineR::PCA(X = Angle_inds,scale.unit = TRUE, 
                          graph = FALSE)
mu_angle<-colMeans(Angle_inds)
percent_variance<-PCA_inds$eig[c(1:2),2]
percent_variance
plot(c(100,PCA_inds$eig[,2]),type="o"
     ,main=paste0("Evolution of the ratio of unexplained inertia for ",ifelse(Nom_v=="Surcote", 
                                                                              "S",Nom_v)))
sd_angle<-apply(X =Angle_inds,MARGIN = 2,FUN = sd)
resultat<-scale(Angle_simul,center = mu_angle,
                scale = sd_angle)
coords_donnees<-PCA_inds$ind$coord[,c(1:2)]
conversion_donnees<-resultat%*%PCA_inds$svd$V[,c(1:2)]
sum(PCA_inds$svd$V[,1]^2)
ks.test(conversion_donnees[,1],coords_donnees[,1])
wilcox.test(conversion_donnees[,1],coords_donnees[,1])


# Keep null hypothesis for the first dimension (KS)---------------------------
# Same thing with Mann Whintey.
ks.test(conversion_donnees[,2],coords_donnees[,2])
wilcox.test(conversion_donnees[,2],coords_donnees[,2])

Indices_for_illustration<-sample(x = c(1:nrow(conversion_donnees)),
                                 size = nrow(coords_donnees))
df_simul_forme<-cbind.data.frame(conversion_donnees)
colnames(df_simul_forme)<-sapply(c(1:ncol(df_simul_forme)), 
                               function(x){return(paste0("Score_",x))})
summary(df_simul_forme)

df_data_forme<-as.data.frame(coords_donnees)
colnames(df_data_forme)<-sapply(c(1:ncol(df_data_forme)), 
                              function(x){return(paste0("Score_",x))})
summary(df_data_forme)
# define a new df to use a categorical variable ---------------------------

df_combo<-rbind.data.frame(df_data_forme,df_simul_forme)
df_combo$Legend<-c(rep("data",nrow(df_data_forme)),
                 rep("simulations",nrow(df_simul_forme)))
df_combo
summary(df_data_forme)

plot(density(df_simul_forme$Score_2))
head(df_combo)
GG1<-ggplot(data=df_combo,aes(x=Score_1,y=Score_2,colour=Legend,
                              shape=Legend))+
  geom_point(aes(size=Legend))+
  scale_size_manual(values=c("simulations"=0.75,"data"=1.5))+
  geom_xsidedensity(data=df_combo,aes(fill=Legend), alpha = 0.5)+
  geom_ysidedensity(data=df_combo,aes(fill=Legend), alpha = 0.5)+
  xlab(paste0("First dimension (",round(percent_variance[1],1),"% of the variance)"))+
  ylab(paste0("Second dimension (",round(percent_variance[2],1),"% of the variance)"))+
  scale_color_manual(values=cols_)+
  scale_fill_manual(values=cols_)+
  scale_shape_manual(values = c("simulations"=17,"data"=19))+
  guides(fill="none")+
  theme(axis.title=element_text(size=15),
        legend.text=element_text(size=10))
GG1
  

# Analyse des dates obtenues.  --------------------------------------------
Dates_exts<-resultat_Simul_Obs$dates_evenements_extremes
Dates_prises<-read.csv(file="../ss_tend/HIVER/dates_prises.csv")[,2]

d_POIXCT<-as.POSIXct(Dates_prises, format="%d/%m/%Y")
Nb_annees<-diff(range(lubridate::year(d_POIXCT)))
NPY<-length(Dates_prises)/Nb_annees

Annee_ext<-substr(Dates_exts,7,10)

# Chargement limites prévues par EV pour les résidus ----------------------
DF_EV_resid<-read.csv(file=paste0("../Travail_hypothèse_RV/EVA_",Nom_v,"_residus.csv"))[,2:5]

Melteur_resid<-melt(t(DF_EV_resid))
colnames(Melteur_resid)<-c("parametre","temps","estimateur")
# Chargement des limites prévues par EV -----------------------------------

Titre_EV_simul<-paste0("Simulations/observations obtenues pour la variable ",Nom_v," avec limites par (EVA)")
Melteur<-melt(t(resultat_Simul_Obs$simulations_AR))
colnames(Melteur)<-c("Temps","Individu","valeur")

M2<-melt(t(resultat_Simul_Obs$observations_exts_base))
colnames(M2)<-c("Temps","Individu","valeur")
Titre_inds<-paste0("Observations extrêmes pour la variable ",Nom_v," avec limites par (EVA)")
min_y<-min(min(M2$valeur),min(Melteur$valeur))
max_y<-max(max(M2$valeur),max(Melteur$valeur))

Melteur$Individu<-as.character(Melteur$Individu)
M2$Individu<-as.character(M2$Individu)
Melteur$Temps<-(Melteur$Temps-19)/6
GG1<-ggplot(data=Melteur,aes(x=Temps,y=valeur,group=interaction(Individu),col=Individu))+
  geom_line()+
  xlab(" ")+ylab(" ")+
  labs(col="Simulation")+
  ylim(c(min_y,max_y))+
  #ggtitle(paste0("Simulations d'individus extrêmes pour ",Nom_v," (n=",nrow(Simul),")"))+
  guides(col="none")
M2$Temps<-(M2$Temps-19)/6
GG2<-ggplot(data=M2,aes(x=Temps,y=valeur,group=interaction(Individu),col=Individu))+
  geom_line()+
  labs(col="Indice_individu")+
  xlab(" ")+ylab(" ")+
  ylim(c(min_y,max_y))+
  #ggtitle(paste0("Observations extrêmes (n=",nrow(Obs_ext),")"))+
  guides(col="none")
yleft <- textGrob(paste0(nom_present," (m)"),
                  rot=90,
                  gp = gpar(col = "black", fontsize = 15))
Xbottom<-textGrob("Time (hour) with respect to tidal peak",
                  gp = gpar(col = "black", fontsize = 15))
GRID_simul_vs_obs<-gridExtra::grid.arrange(GG1,GG2,ncol=2,
                        bottom=Xbottom,
                        left=yleft)
ggsave(GRID_simul_vs_obs,
       filename="graphiques_Surcote/exs_simuls_vs_obs.png",
       width = 12,
       height=7)

# Visualisation_Johanna_vs_preds ------------------------------------------
#########
Predictions_for_Johanna<-resultat_Simul_Obs$liste_Pred_Johanna$simulations
courbe_Johanna<-resultat_Simul_Obs$liste_Pred_Johanna$X1
Origine<-resultat_Simul_Obs$liste_Pred_Johanna$X0

# (a) Comparing_values ------------------------------------------------------
Courbes_predites<-melt(t(Predictions_for_Johanna))
colnames(Courbes_predites)<-c("Time","ind","value")
Courbes_predites$ind<-as.character(Courbes_predites$ind)
Courbes_predites$realite_Joh<-as.numeric(courbe_Johanna)[Courbes_predites$Time]
Courbes_predites$origine<-as.numeric(Origine)[Courbes_predites$Time]

ggplot(data=Courbes_predites,aes(x=Time,y=value,group=interaction(ind),
                                 col=ind,linetype="simulations"))+
  geom_line(alpha=0.3)+
  guides(col="none")+
  geom_line(aes(x=Time,y=realite_Joh,linetype="Johanna observation"))+
  geom_line(aes(x=Time,y=origine,linetype="previous_observation"))+
  labs(linetype="Legend")+
  ylab("value (m)")+
  ggtitle(paste0("Predictions for Johanna observation of the ",
                 resultat_Simul_Obs$liste_Pred_Johanna$date_Johanna," for ",
                 ifelse(Nom_v=="Surcote","S",Nom_v))) +
  annotate("text", x = 30, y = as.numeric(Origine[30]), 
           label=resultat_Simul_Obs$liste_Pred_Johanna$date_X0)

nom_recherche<-paste0("../ss_tend/HIVER/Surcote_ss_tend.csv")
obs<-read.csv(nom_recherche)[,2:38]
colnames(obs)<-c(1:37)
Dens<-melt(t(obs))
colnames(Dens)<-c("time","ind","value")
Dens$time<-as.factor(Dens$time)
ggplot(data=Dens,aes(x=time,y=value,fill=time))+
  geom_boxplot(width=0.5)+
  geom_violin(show.legend = FALSE)
  
Inds_exts<-resultat_Simul_Obs$indices_exts
Dens2<-melt(t(obs[Inds_exts,]))
colnames(Dens2)<-c("time","ind","value")
Dens2$time<-(Dens2$time-19)/6
Dens2$time<-as.factor(Dens2$time)
ggplot(data=Dens2,aes(x=time,y=value,fill=time))+
  geom_boxplot(width=0.5)+
  geom_violin(show.legend = FALSE)+
  ylab("Surge (m)")
summary(Dens2)
# (b) Comparing the shape -------------------------------------------------

# First version --------------------------------------------------------
Pred_Joh_std<-scale(Predictions_for_Johanna,center=mu,
                    scale=sd)
Coords_pred_Joh_brut<-Pred_Joh_std%*%PCA_brut_$svd$V[,c(1:2)]

Joh_std<-((as.numeric(courbe_Johanna))-mu)/sd
Coords_Joh_brut<-Joh_std%*%PCA_brut_$svd$V[,c(1:2)]

plot(Coords_pred_Joh_brut,main="Coordinated predictions vs reality (cross in red) - - ACP")
points(Coords_Joh_brut,col="red",pch=4,cex=2)

# sweep to remove a vector from a matrix ----------------------------------

# 2 to give what explains at each column ----------------------------------

Distance_each_coord<-sweep(Coords_pred_Joh_brut, 2, 
                               Coords_Joh_brut)^2
Distance_to_reality<-apply(Distance_each_coord,MARGIN = 1,FUN = sum)
summary(Distance_to_reality)

# comparing with general  -------------------------------------------------
Distance_each_coord_for_allS<-sweep(resultat_brut%*%PCA_brut_$svd$V[,c(1:2)], 2, 
                                    Coords_Joh_brut)^2
Distance_to_reality_allS<-apply(Distance_each_coord_for_allS,MARGIN = 1,
                                FUN = sum)
summary(Distance_to_reality_allS)

NB_knn<-10
indices_d_min<-order(Distance_to_reality)[1:NB_knn]
points(Coords_pred_Joh_brut[indices_d_min,],col="orange")
matplot(t(Predictions_for_Johanna[indices_d_min,]), 
        type="l")
lines(x=c(1:37),Joh)
df_sub<-melt(t(Predictions_for_Johanna[indices_d_min,]))
colnames(df_sub)<-c("Time","ind","value")
df_sub$ind<-as.character(df_sub$ind)
df_sub$Joh<-Joh[df_sub$Time]
ggplot(data=df_sub,aes(x=Time,y=value,col=ind,group=interaction(ind),linetype="nearest_simul"))+
  geom_line()+
  geom_line(aes(x=Time,y=Joh,linetype="Johanna"))+
  guides(col="none")+
  labs(linetype="Legend")+
  ggtitle(paste0("Nearest simulated time series to Johanna's for ",Nom_v," (K=",NB_knn,")"))

# variante_distance_L2 -------------------------------------------------------
Distance_by_time<-sweep(Predictions_for_Johanna,2,
      Joh)
Distance_L2<-apply(abs(Distance_by_time),MARGIN = 1,FUN = calcul_norme_L2)

indices_d_minL2<-order(Distance_L2)[1:NB_knn]
matplot(t(Predictions_for_Johanna[indices_d_minL2,]),type="l")
lines(x=c(1:37),Joh)
print(indices_d_min)
print(indices_d_minL2)
print("Même indices norme L2 coords et L2")

# geom_point(aes(x =Coords_Joh_brut[,1], y = Coords_Joh_brut[,2], color = "Annotated Point"), size = 3)
# ML on reconstructed trajectories ---------------------------------
### Resul ML -------------------------------------------------------------")
N_angle<-paste0("Angle ",Nom_v)
source("../fonctions/fonctions_perfs_ML.R")
racine_ML<-c("resultat_ML/residus/")

# SVM ---------------------------------------------------------------------
ENSEMBLE_liens<-c("linear","sigmoid","polynomial","radial")
ech_train<-list("simul"=Simulations[Indices_train_s,],
                "real"=Realite[Indices_train_r,])
ech_test<-list("simul"=Simulations[Indices_test_s,],
               "real"=Realite[Indices_test_r,])

# Graphiques_pr_test ------------------------------------------------------

Melteur<-melt(t(Simulations[Indices_test_s,]))
colnames(Melteur)<-c("Temps","Individu","valeur")
Melteur$Individu<-as.character(Melteur$Individu)

Q1_reality<-apply(X = Realite,MARGIN = 2,FUN = function(x){
  return(as.numeric(quantile(x,0.05)))
})
Q2_reality<-apply(X = Realite,MARGIN = 2,FUN = function(x){
  return(as.numeric(quantile(x,0.95)))
})

Q1_simul<-apply(X = Simul,MARGIN = 2,FUN = function(x){
  return(as.numeric(quantile(x,0.05)))
})
Q2_simul<-apply(X = Simul,MARGIN = 2,FUN = function(x){
  return(as.numeric(quantile(x,0.95)))
})


M2<-melt(t(Realite[Indices_test_r,]))
colnames(M2)<-c("Temps","Individu","valeur")
M2$Individu<-as.character(M2$Individu)

M2$Q05<-Q1_reality[M2$Temps]
M2$Q95<-Q2_reality[M2$Temps]
Melteur$Q05<-Q1_simul[Melteur$Temps]
Melteur$Q95<-Q2_simul[Melteur$Temps]
Titre_inds<-paste0("Observations extrêmes pour la variable ",Nom_v," avec limites par (EVA)")
min_y<-min(min(M2$valeur),min(Melteur$valeur))
max_y<-max(max(M2$valeur),max(Melteur$valeur))
Melteur$mu<-mu1[Melteur$Temps]
M2$mu<-mu2[M2$Temps]
M2$Temps<-(M2$Temps-19)/6
Melteur$Temps<-(Melteur$Temps-19)/6

GG1<-ggplot(data=Melteur,
            aes(x=Temps,y=valeur,col=Individu,
                             group=interaction(Individu), 
                             ))+
  geom_ribbon(aes(ymin=Q05,ymax=Q95,linetype="confidence_band"),
             alpha=0.02,fill="grey",col="darkblue")+
  guides(col="none")+
  geom_line(aes(x = Temps,y=mu,linetype="mean"),
            col="red",size=1.1)+
  geom_line(alpha=0.7)+
  scale_linetype_manual("Legend",values=c("confidence_band"=2,
                                            "mean"=5))+
  theme(axis.title=element_text(size=15),
        legend.text=element_text(size=10))+
  xlab(" ")+
  ylab(" ")+
  ylim(c(min_y,max_y))

GG1

GG2<-ggplot(data=M2,aes(x=Temps,y=valeur,col=Individu,group=interaction(Individu)))+
  geom_ribbon(mapping=aes(ymin=Q05,ymax=Q95,linetype="confidence_band"),
              alpha=0.02,
              fill="grey",col="darkblue")+ 
  geom_line(aes(y=mu,linetype="mean"),col="red",
            size=1.1)+
  geom_line(alpha=0.7)+
  guides(col="none")+
  theme(axis.title=element_text(size=15),
        legend.text=element_text(size=10))+
  scale_linetype_manual("Legend",values=c("confidence_band"=2,
                                        "mean"=5))+
  xlab(" ")+
  ylab(" ")+
  ylim(c(min_y,max_y))
GG2

gridExtra::grid.arrange(GG1,GG2,ncol=2, 
                        bottom=Xbottom,
                        left=yleft)

LISTE_ML<-fnct_ML_realite_simulations_all_kernels(ech_test = ech_test,
                                              ech_train = ech_train,
                                                   variable = Nom_v)

# Pour l'angle ------------------------------------------------------------
Cas_non_nul_train<-intersect(inds_non_nuls,Indices_train_s)
Cas_non_nul_test<-intersect(inds_non_nuls,Indices_test_s)

Train_simul_angle<-t(t(Simul[Cas_non_nul_train,])%*%diag(NL2_simul[Cas_non_nul_train]^(-1)))
Test_simul_angle<-t(t(Simul[Cas_non_nul_test,])%*%diag(NL2_simul[Cas_non_nul_test]^(-1)))
ech_train<-list("simul"=Train_simul_angle,
                "real"=Angle_inds[Indices_train_r,])
ech_test<-list("simul"=Test_simul_angle,
               "real"=Angle_inds[Indices_test_r,])
N_angle<-paste0("Angle de ",Nom_v)

LISTE_ML_angle<-fnct_ML_realite_simulations_all_kernels(ech_test = ech_test,
                                                  ech_train = ech_train,
                                                  variable = N_angle)


# Perfs with several draws -----------------------------------------------------------
########################

H_Params<-c("radial","logit",50)
NB_times<-100
ALPHA_PROP<-0.10
K<-4
typeS<-"simple"
resultat<-matrix(NA,ncol=3,nrow=3)

# Use of the function -----------------------------------------------------
Lperfs_Y<-Running_perfs_ML(Base_simul = Simul,Base_data = Realite, 
                   hyp_param = H_Params, NB_times =100,alpha_prop=ALPHA_PROP,K = K, 
                   type_sampling = typeS,title_ROC =paste0(nom_present),
                   NB_shown_ROC = 20)

Lperfs_Y$HSS
resultat[1,]<-Lperfs_Y$Accuracy$qu_accuracy

Lperfs_L2<-Running_perfs_ML(Base_simul = NL2_simul,Base_data = NL2_ext, 
                              hyp_param = H_Params, NB_times = NB_times,alpha_prop=ALPHA_PROP, 
                              K=K,type_sampling = typeS
                           ,title_ROC = paste0(nom_present," Cost-f"),NB_shown_ROC = 20)

Lperfs_L2$HSS
resultat[2,]<-Lperfs_L2$Accuracy$qu_accuracy

Lperfs_A_f<-Running_perfs_ML(Base_simul = Angle_simul,Base_data = Angle_inds, 
                              hyp_param =H_Params , NB_times = NB_times,alpha_prop=ALPHA_PROP, 
                              K=K,type_sampling = typeS, 
                            title_ROC = paste0(nom_present," angle"),NB_shown_ROC = 20)
Lperfs_A_f$HSS
resultat[3,]<-Lperfs_A_f$Accuracy$qu_accuracy
resultat

# Export des resultats ------------------------------------------------------------------
resultat
df<-rbind.data.frame(H_Params,resultat)
df[1,]<-paste0(c("kernel : ","link : ","Nb-trees : "),df[1,])
rownames(df)<-c("Hyper-params","X","l(X)","A(X)")
colnames(df)<-c(paste0("SVM"),paste0("GLM"),"RForest")
df<-tibble::rownames_to_column(df, "Input")
df
write.csv(x = df,file = paste0("resultat_ML/residus/resultat_ML_",Nom_v,"_",methode_,"_",typeS,".csv"), 
          row.names=FALSE,quote = FALSE)


# L2 ----------------------------------------------------------------------
ech_train<-list("simul"=NL2_simul[Indices_train_s],
                "real"=NL2_ext[Indices_train_r])
ech_test<-list("simul"=NL2_simul[Indices_test_s],
               "real"= NL2_ext[Indices_test_r])
LISTE_ML_L2<-fnct_ML_realite_simulations_all_kernels(ech_test = ech_test,
                                                        ech_train = ech_train,
                                                        variable = paste0("L2 de ",Nom_v))
l_kernel<-names(LISTE_ML)[c(1:4)]
for(kernel in l_kernel){
  x_vect<-rbind(LISTE_ML[[kernel]]$resultat_conf_test,
                LISTE_ML_L2[[kernel]]$resultat_conf_test,
                LISTE_ML_angle[[kernel]]$resultat_conf_test)
                
  if(!is.null(x_vect)){
    write.csv(paste0("resultat_ML/residus/svm_individu/",kernel,"_",Nom_v,"_",methode_,".csv"),x=x_vect)
  }
}

# GLM avec fonctions de lien ----------------------------------------------
ech_train<-list("simul"=Simulations[Indices_train_s,],
                "real"=Realite[Indices_train_r,])
ech_test<-list("simul"=Simulations[Indices_test_s,],
               "real"=Realite[Indices_test_r,])
Liste_fonction_lien<-c("logit", "probit", "cloglog","log","cauchit")
LISTE_resultat_glm<-fnct_logistigue_glm_mult_links(ech_test = ech_test,
                                                   ech_train = ech_train,
                                                   liste_fncts = Liste_fonction_lien,
                                                   variable = Nom_v)




# Travail angle -----------------------------------------------------------

ech_train<-list("simul"=Train_simul_angle,
                "real"=Angle_inds[Indices_train_r,])
ech_test<-list("simul"=Test_simul_angle,
               "real"=Angle_inds[Indices_test_r,])
LISTE_angle_glm<-fnct_logistigue_glm_mult_links(ech_test = ech_test,
                                                ech_train = ech_train,
                                                liste_fncts = Liste_fonction_lien,
                                                variable = N_angle)
matplot(t(Train_simul_angle),type="l")
matplot(t(Test_simul_angle),type="l")

# Travail L2 --------------------------------------------------------------

ech_train<-list("simul"=NL2_simul[Indices_train_s],
                "real"=NL2_ext[Indices_train_r])
ech_test<-list("simul"=NL2_simul[Indices_test_s],
               "real"= NL2_ext[Indices_test_r])

LISTE_glm_L2<-fnct_logistigue_glm_mult_links(ech_test = ech_test,
                                         ech_train = ech_train,
                                         liste_fncts = Liste_fonction_lien,
                                         variable=paste0("L2 de ",Nom_v))
for(lien in Liste_fonction_lien){
  x_vect<-rbind(LISTE_resultat_glm[[lien]]$resultat_conf_test,
                LISTE_glm_L2[[lien]]$resultat_conf_test,
                LISTE_angle_glm[[lien]]$resultat_conf_test
                )
  if(!is.null(x_vect)){
    write.csv(paste0("resultat_ML/residus/glm_individu/",lien,"_",Nom_v,"_",methode_,".csv"),x=x_vect)
  }

}

# Export des performances des modèles -------------------------------------
##############

svm_perf<-LISTE_ML$performances_globales
glm_perf<-LISTE_resultat_glm$performances_globales
write.csv(paste0("resultat_ML/residus/glm_individu/performances_glm_",Nom_v,".csv"),x=glm_perf)
write.csv(paste0("resultat_ML/residus/svm_individu/performances_svm_",Nom_v,".csv"),x=svm_perf)

# Les résultats dépendent de ce qu'on compare -----------------------------
# Angle versus les individus.


# Analyse de la loi de la norme -------------------------------------------
##################

norme_<-c(NL2_simul,NL2_ext)
origine<-c(rep("simulations",length(NL2_simul)),rep("data",length(NL2_ext)))
df<-cbind.data.frame(norme_,origine)
colnames(df)<-c("Norm","Origin")

GG_norme<-ggplot(data=df,aes(y=Norm,x=Origin,fill=Origin))+
  geom_violin(alpha=0.6)+
  geom_boxplot(width=0.3)+
  ylab("L2 norm (m)")+
  #ggtitle(paste0("Distribution of the L2 norm of ",Nom_graph))+
  labs(fill="Legend")+
  scale_fill_manual(values = cols_)+
  theme(axis.title=element_text(size=15),
        legend.text=element_text(size=10))
print(GG_norme)

# Niveau de retour --------------------------------------------------------
NL2<-apply(resultat_Simul_Obs$X,MARGIN = 1,FUN = calcul_norme_L2)
inds_exts_ou_non<-(c(1:length(NL2))%in%Indices_exts)
seuil<-max(NL2[which(NL2<min(NL2_ext))])

NL2_ext_tf<-NL2[which(NL2>=min(NL2_ext))]
Modele_ev_dummy<-extRemes::fevd(x = NL2,threshold = seuil,type="GP")
rate_exceedance<-length(NL2_ext_tf)/length(NL2)

# ss_ech que l'on a ... ---------------------------------------------------
goftest::ad.test(x = NL2_ext_tf,null = extRemes::"pevd",scale=Modele_ev_dummy$results$par[1],
                 shape=Modele_ev_dummy$results$par[2],threshold=seuil,
                 type="GP")

periods_years<-c(2,5,10,20,50,80,
                 100,120,200,
                 250,300,500,
                 800)
probabilite<-1-(1/periods_years*rate_exceedance)
q_simulations<-quantile(NL2_simul,probabilite)
estimateur_empirique<-as.numeric(quantile(NL2_ext,probabilite))

# Boostrap pour identifier les limites ------------------------------------------------
Fnct<-function(x,taille_ech,proba_retour){
  Ech<-sample(c(1:length(x)),size = 100,replace = TRUE)
  x_tilde<-x[Ech]
  q_estim<-as.numeric(quantile(x_tilde,proba_retour))
  return(q_estim)
}
Fnct_par_periode<-function(series,taille_ech,periode,B){
  Echantillons_une_periode<-replicate(n =B,
                                      expr = Fnct(x=series,
                                                  taille_ech =taille_ech ,
                                                  proba_retour=periode))
  estimateur<-as.numeric(quantile(series,periode))
  Intervalles_1_periode<-as.numeric(quantile(Echantillons_une_periode,probs = c(0.025,0.975)))
  return(c(Intervalles_1_periode,estimateur))
}

# Import date correspondante ----------------------------------------------

shape<-Modele_ev_dummy$results$par[[2]]
scale<-Modele_ev_dummy$results$par[[1]]

estimateur_EV<-extRemes::rlevd(period = periods_years,threshold  = seuil,
                scale = scale,
                shape=shape, 
                type="GP",
                rate = rate_exceedance,
                npy=NPY)
periods_years_extended<-c(seq.int(from = 0.2,min(periods_years),by=0.2),periods_years)
periods_years_extended
m<-(periods_years_extended*NPY)
Prop<-m*rate_exceedance
entree<-1-(1/Prop)
Niveaux_empiriques_simul<-as.numeric(quantile(x = NL2_simul,entree))
simul_<-cbind.data.frame(periods_years_extended,Niveaux_empiriques_simul)
colnames(simul_)<-c("p_years","simul_niveau")

RL_plot<-plot(Modele_ev_dummy,type = "rl",
     main=paste0("Niveau de retour par EVA pour la norme ",Nom_v),
     rperiods=c(periods_years))
empirique<-RL_plot$empirical

Modele_vs_emp<-RL_plot$model

# Johanna jour même (vs) simulations -------------------------------------------------------
Simuls_obtenues<-resultat_Simul_Obs$liste_Johanna$simulations
Vraie_obs<-resultat_Simul_Obs$liste_Johanna$vraie_obs
Q2<-apply(X = Simuls_obtenues,MARGIN = 2,FUN = function(x){return(quantile(x,0.975))})
Q1<-apply(X = Simuls_obtenues,MARGIN = 2,FUN = function(x){return(quantile(x,0.025))})
Max_g<-max(max(Q1),max(Q2),max(Vraie_obs))
Min_g<-min(min(Q1),min(Q2),min(Vraie_obs))


# analyse EV --------------------------------------------------------------

simul<-resultat_Simul_Obs$simulations_AR

periods_years<-c(2,5,10,20,50,80,
                 100,120,200,
                 250,300,500,
                 800)
NB_annee<-37
limite_sup<-3*37
limite_sup
periods_years<-periods_years[which(periods_years<limite_sup)]
periods_years

Dates<-read.csv(file="../ss_tend/HIVER/dates_prises.csv")[,2]

cols_<-c("data"="blue","simulations"="orange","confidence_band"="darkblue")
p_u<-rep(0.05,ncol(simul))
fich_inds<-read.csv(paste0("residus_clust/",inds_extremes))
inds_exts<-fich_inds$Inds_Exts
alpha<-0.05
if(Nom_v=="Hs"){
  Y_graph_lim<-c(4,11)
}
if(Nom_v=="Surcote"){
  Y_graph_lim<-c(0.13,1)
}

for(t in c(13,19,25,31)){
  lag<-t-19
  comment<-ifelse(lag<0,paste0(abs(round(lag/6)),"h before"),
                  ifelse(lag==0," at",paste0(round(lag/6),"h after")))
  series_simul<-simul[,t]
  series_t<-resultat_Simul_Obs$X[,t]
  seuil<-quantile(series_t,1-p_u[t])
  
  series_ext<-subset(series_t,series_t>seuil)
  NB_years_total<-37
  npy<-length(series_t)/NB_years_total
  p_proba<-1-(npy*periods_years*p_u[t])^(-1)
  Nom_graph<-ifelse(nom_present=="Surcote","S",
                    nom_present)
  # origine : package extRemes.  --------------------------------------------
  RL_ggplot_cond_ext(series = series_t,seuil = seuil,
                     period_years = periods_years,
                     NPY = round(npy),
                     plus_simul = TRUE,
                     series_simul = series_simul, 
                     titre = paste0("Return level (ML) ",comment, " the tidal peak for ",Nom_graph), 
                     nom_variable = Nom_v, 
                     cols_ggplot = cols_, 
                     unit_used = "m",
                     ylim_opt = Y_graph_lim, 
                     alpha=0.05,
                     Individus_exts=Indices_exts)
  # Niveau retour pr les donnees --------------------------------------------
  xp2<-ppoints(length(series_t))
  ysorted<-sort(series_t)
  tf_period_data<--1/log(xp2)[ysorted> seuil]/npy
}
summary(series_simul)


# Fin du code -------------------------------------------------------------
# Graphiques --------------------------------------------------------------


