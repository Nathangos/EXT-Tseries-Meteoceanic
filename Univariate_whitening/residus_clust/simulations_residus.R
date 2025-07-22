rm(list=ls())
set.seed(133)
par(mfrow=c(1,1))

#### Import des packages. #####
source("../fonctions/fonctions_Pareto.R")
source("../fonctions/fonctions.R")
require(dplyr)
require(ggplot2)
require(reshape2)
require(FactoMineR)
require(extRemes)
require(fda)
require(parallel)
require(VineCopula)
require(grid)
require(gridExtra)
require(mvPot)

### Source des fonctions. ####
F_names<-ls()
nb_coeurs<-detectCores()-4


# Ouverture des clusters --------------------------------------------------
coeurs<-makeCluster(nb_coeurs)
clusterExport(coeurs,varlist=F_names)

Simulation_ACP_TS_EXT<-function(type_donnees,coeurs,NB_dim_ACP,lien_donnees,nom_variable,p_U,p_L2,n.dens,opt_Frech,M,std_donnees,type_entree,type_vario){
  if(type_entree=="clust"){
    All<-read.csv(file=lien_donnees)[,2:39]
    colnames(All)<-c(1:38)
    Vrais_indices<-All[,38]
    Donnes<-All[,c(1:37)]
  }
  else{
    Donnes<-read.csv(file=lien_donnees)[,2:38]
  }
  colnames(Donnes)<-c(1:37)
  rownames(Donnes)<-c(1:nrow(Donnes))
  dates_import<-read.csv(file="../ss_tend/HIVER/dates_prises.csv")[,2]
  d_POIXCT<-as.POSIXct(dates_import, format="%d/%m/%Y")
  Nb_annees<-diff(range(lubridate::year(d_POIXCT)))
  NPY<-nrow(Donnes)/Nb_annees
  AN_GPD<-Analyse_seuil_GPD(donnees = Donnes,fonction_seuil = p_U,n.dens = n.dens,
                            nom=nom_variable,type_entree=type_entree,
                            dates_prises=dates_import,
                            j_show=c(19))
  if(type_entree=="donnees_brutes"){
    return(TRUE)
  }
  Vecteur_DIM<-c(1:ncol(Donnes))

  Resultat_P<-as.data.frame(t(sapply(Vecteur_DIM,f_marginales_all_Pareto,data_to_tf = Donnes,p_u=p_U,
                                     n.dens=n.dens)))
  K<-Resultat_P$kernel_dens
  PARETO<-Resultat_P$obs
  L_EVT<-list("gamma"=Resultat_P$gamma,
              "scale"=Resultat_P$scale,
              "threshold"=Resultat_P$threshold,
              "p_u"=Resultat_P$p_u)
  vecteur_pareto<-do.call(cbind.data.frame,PARETO)

  # Hypothèse Pareto pour la transformation ---------------------------------
  Analyse_Pareto_par_temps(donnees_Pareto = vecteur_pareto,nom = nom_variable,n.dens=n.dens)
  
  if(opt_Frech==TRUE){
    Unif<-1-(as.matrix(vecteur_pareto))^(-1)
    #Ici frechet, inverse fonction de repartition d'une Frechet.
    vecteur_pareto<-(-1)*(log(Unif))^(-1)
    for(t in c(1:ncol(vecteur_pareto))){
      serie_frech_t<-vecteur_pareto[,t]
      #POT::tsdep.plot(serie_frech_t,u = 0.95,lag.max = 10,
      #               main=paste0("Resultat pour le temps pour T(",nom_variable,") en t=",t))
      # Picklands ---------------------------------------------------------------
      h<-c(1:50)
      # fnct_Picklands_lag_h<-function(h,series){
      #   fin<-length(series)-h
      #   x_<-series[c(1:fin)]
      #   debut<-h+1
      #   x_lag1<-series[c(debut:length(series))]
      #   u1<-quantile(x_,0.90)
      #   u2<-quantile(x_lag1,0.90)
      #   modele<-POT::fitbvgpd(data = cbind(x_,x_lag1),
      #                 threshold = c(u1,u2), 
      #                 model = "log",
      #                 alpha=1)
      #   POT::pickdep(object = modele,main=paste0("Pickdep at t=",t," h=",h))
      #   
      # }
      # sapply(X = h,FUN = fnct_Picklands_lag_h,series=serie_frech_t)
    }
  }
  Vecteur_NL2_original<-apply(X = Donnes,MARGIN = 1,
                                 FUN = calcul_norme_L2)
  
  LNL2<-600
  vect_k<-c(20:LNL2)

  # Intervalle de confiance -------------------------------------------------
  # Intervalles_conf<-sapply(X = vect_k,FUN = function_ML_extRemes,NB_years=Nb_annees,
  #                          data_d=Vecteur_NL2_original)
  # plot_ML_int<-fonction_MLplot_resume(resultatML = Intervalles_conf,vecteur_k = vect_k,
  #                                     nom_variable = paste0("Norme L2 de ",nom_variable),
  #                        lims_Y = c(-1,1))
  # print(plot_ML_int)
  
  # Avec autres estimateurs -------------------------------------------------
  # Graphics_estimators_gamma(series = Vecteur_NL2_original,vect_k = vect_k,
  #                              Title_graphic = paste0("Shape parameter of ",nom_variable," L2 norm"),
  #                              NB_years=Nb_annees)
  
  # Analyse GPD norme ech orig ----------------------------------------------
  
  # Vérifier que les excédents à l'echelle  ---------------------------------
  # suivent bien une loi GPD. -----------------------------------------------

  seuil<-as.numeric(quantile(Vecteur_NL2_original,1-p_L2))
  # Outils_POT_graphique(seuil = seuil,Q1 = 0.50,
  #                      Q2 = 0.98,series = Vecteur_NL2_original,
  #                      dates = dates_import,
  #                      titre_variable=paste0("Norme L2 de ",nom_variable))
  modele_ev_orig<-fevd(x = Vecteur_NL2_original,threshold = seuil,
       type = "GP")
  
  vecteur_norm_L2<-apply(X= vecteur_pareto,MARGIN = 1,FUN =calcul_norme_L2)
  # Graphics_estimators_gamma(series = vecteur_norm_L2,vect_k = vect_k,
  #                              Title_graphic = paste0("Shape parameter of T(",nom_variable,")'s L2 norm"),
  #                              NB_years = Nb_annees)
  
  rownames(vecteur_pareto)<-1:nrow(vecteur_pareto)
  colnames(vecteur_pareto)<-c(1:ncol(vecteur_pareto))
  
  v_seuil<-quantile(vecteur_norm_L2,1-p_L2)
  #v_seuil<-40
  print(mean(as.numeric(vecteur_norm_L2>v_seuil)))
  
  Inds_Exts<-which(vecteur_norm_L2>v_seuil)
  Inds_non_exts<-which(vecteur_norm_L2<=v_seuil)
  matplot(t(Donnes[Inds_Exts,]),type="l")
  
  seuil_transfere<-min(Vecteur_NL2_original[Inds_Exts])
  mtext(paste0("Observations extrêmes conservées(",nom_variable,")"))
  Exts<-vecteur_pareto[Inds_Exts,]
  
  # Observation des extrêmes pour les scores --------------------------------
  # une fois la conversion en M-S simple ------------------------------------
  L<-nrow(vecteur_pareto)/10
  Nom_alt<-paste0("T(",nom_variable,")")
  Analyse_extreme_proj(base_RV = vecteur_pareto,nb_scores = NB_dim_ACP,
                        stand_ = std_donnees,nom_variable = Nom_alt,
                        L=L,NB_year = Nb_annees,M1=200,M2=300)

  # Composantes angulaires --------------------------------------------------
  Vect_L2_exts<-vecteur_norm_L2[Inds_Exts]
  KS_norme<-ks.test(x = Vect_L2_exts,extRemes::"pevd",
          threshold=v_seuil,
          scale=v_seuil,
          shape=1,type="GP")
  
  AD_norme<-goftest::ad.test(x = Vect_L2_exts,extRemes::"pevd",
                    threshold=v_seuil,
                    scale=v_seuil,
                    shape=1,type="GP")
  print("###### Test statistique KS/Anderson Darling (norme L2/seuil) suit une Pareto #####")
  print(c(KS_norme$p.value,AD_norme$p.value))
  
  Comp_angulaire<-t(t(Exts)%*%diag(Vect_L2_exts^(-1)))
  matplot(t(Comp_angulaire),type="l")
  mtext(paste0("Composantes angulaires (",nom_variable,")"))

  # convert to log to use after exponential ---------------------------------
  
  Comp_angulaire<-Comp_angulaire
  # ACP --------------------------------------------------
  sd_temps<-as.numeric(apply(X = Comp_angulaire,
                             MARGIN = 2,FUN = sd))
  esp_temps<-as.numeric(colMeans(Comp_angulaire))
  print("moy par temps")
  print(esp_temps)
  Analyse<-FactoMineR::PCA(X =Comp_angulaire,
                           scale.unit = std_donnees,
                           ncp = NB_dim_ACP)
  val_lambda<-Analyse$eig[,1]
  vect_diff<-cumsum(c(0,val_lambda))/sum(val_lambda)
  vecteur_propvarexp<-1-vect_diff
  L<-length(vecteur_propvarexp)-1
  plot(c(0:L),vecteur_propvarexp,
       type="o",ylab="Proportion of unexplained inertia", 
       xlab="Number of eigenvectors",
       cex.lab=1.5)
  par(mfrow=c(1,2))
  Coord1<-Analyse$ind$coord[,1]
  Coord2<-Analyse$ind$coord[,2]
  F_propres<-Analyse$svd$V
  plot(Coord1,Coord2)
  Coordonnees<-Analyse$ind$coord[,1:NB_dim_ACP]
  write.csv(file=paste0("residus_clust/ACP_results/Theta_coords_ACP_",nom_variable,".csv"),
            x=Coordonnees)
  dataUnif<-pobs(Coordonnees)
  # Simulation_norme --------------------------------------------------------
  valeurs_unif<-runif(M)
  reals_P_std<-sapply(valeurs_unif,Generation_Pareto_std)
  ECHELLE<-reals_P_std*v_seuil
  
  if(NB_dim_ACP==2){
    Famille_choisie<-BiCopSelect(dataUnif[,1],dataUnif[,2])
    Realisations_copules<-BiCopSim(N=M,family =Famille_choisie$family,par = Famille_choisie$par,par2 =Famille_choisie$par2,check.pars = TRUE)
    sommaire<-summary(Famille_choisie)
    
    Short_famille_name<-BiCopName(family = sommaire[["family"]],short=TRUE)
    elt<-c(sommaire[["familyname"]],sommaire[["par"]],sommaire[["par2"]],sommaire[["tau"]],sommaire[["beta"]],unlist(sommaire[["taildep"]]))
    L<-length(elt)
    elt<-as.data.frame(t(elt))
    colnames(elt)<-c("Nom","par","par2","tau","beta","ltd","uptd")
    elt[,c("par","par2","tau","beta","ltd","uptd")]<-round(as.numeric(elt[,c("par","par2","tau","beta","ltd","uptd")]),1)
    write.csv(elt,file = paste0("residus_clust/",nom_variable,"_coords_UV_copules.csv"),
             row.names = FALSE)
  }
  else{
    Matrice_C<-function_Structure_Matrice(NB_dim = NB_dim_ACP)
    Famille_choisie<-VineCopula::RVineCopSelect(dataUnif[,c(1:NB_dim_ACP)],
                                                Matrix = Matrice_C)
    Code_utilisee<-summary(Famille_choisie)$family
    write.csv(x = Code_utilisee,
              file = paste0("residus_clust/ACP_results/",nom_variable,"_modele_Vine.txt"))
    Realisations_copules<-VineCopula::RVineSim(N=M,RVM = Famille_choisie) 
    obj_famille<-summary(Famille_choisie)
    obj_famille$par<-ifelse(obj_famille$cop=="I",NA,obj_famille$par)
    obj_famille$par2<-ifelse(obj_famille$cop=="I",NA,obj_famille$par2)
    obj_famille$utd<-ifelse(obj_famille$cop %in% c("I","BB8_90","Tawn2_90","BB8","J270","J90","F"),NA,
                            obj_famille$utd)
    obj_famille$ltd<-ifelse(obj_famille$cop %in% c("I","BB8_90","Tawn2_90","BB8","J270","J90","F","Tawn","J","Tawn2","SC"),NA,
                            obj_famille$ltd)
    obj_famille$family<-BiCopName(obj_famille$family,short = FALSE)
    dim_variable<-ncol(obj_famille)
    obj_famille<-obj_famille[,c(1:3,5:dim_variable)]
    colnames(obj_famille)[c(1:3)]<-c("arbre","couples-cond","nom-copule")
    obj_famille[,c("par","par2","tau","utd","ltd")]<-round(obj_famille[,c("par","par2","tau","utd","ltd")],1)
    write.csv(obj_famille,file = paste0("residus_clust/",nom_variable,"_coords_UV_copules.csv"))
  }
  par(mfrow=c(1,2))
  plot(dataUnif[,1],dataUnif[,2],xlab="Première coordonnée (observations)",ylab="Seconde coordonnée (observations)")
  plot(Realisations_copules[,1],Realisations_copules[,2],col="red",xlab="Première coordonnée",ylab="Seconde coordonnée")
  title(paste0("Comparaison entre les observations et les simulations pour ", nom_variable," en ",tolower(type_donnees)),outer = TRUE,line=-3)
  
  # Conversion dans la bonne échelle ----------------------------------------
  Matrice_coords_ALL<-matrix(NA,ncol = NB_dim_ACP,nrow = M)
  for(l in c(1:NB_dim_ACP)){
    Matrice_coords_ALL[,l]<-quantile(Coordonnees[,l],probs=Realisations_copules[,l])
  }
  if(std_donnees==TRUE){
    FORME_ACP<-function_reconstitution_trajectory_std(Vector_coords = Matrice_coords_ALL,
                                                      Base_functions_p = F_propres,
                                                      NB_dim =NB_dim_ACP,
                                                      mu_t = esp_temps,
                                                      sd_t = sd_temps)
  }
  else{
    FORME_ACP<-fonction_reconstitution_trajectoire(Vecteur_coords = Matrice_coords_ALL,Base_fonctions_p = F_propres,
                                                   NB_dim = NB_dim_ACP,mu_t = esp_temps)
  }
  # MAJ en mettant à 1 la norme L2 ------------------------------------------
  #FORME_ACP<-exp(FORME_ACP)
  L2_proche_1<-apply(X =FORME_ACP,MARGIN = 1,FUN = calcul_norme_L2)
  FORME_ACP<-t(t(FORME_ACP)%*%diag(L2_proche_1^(-1)))
  
  Nouvelles_obs_P_acp<-t(t(FORME_ACP)%*%diag(ECHELLE))
  if(opt_Frech==TRUE){
    indicatrice_pos<-which(apply(X=Nouvelles_obs_P_acp,
                                 FUN=fonction_trajectoire_positive,MARGIN = 1)==TRUE)
    Unif<-exp(-Nouvelles_obs_P_acp[indicatrice_pos,]^(-1))
    PTO<-(1-Unif)^(-1)
    #seuil_marg<-1+10^(-5)
    seuil_marg<-log(10^5)^(-1)
    indicatrice_pos2<-which(apply(X=(Nouvelles_obs_P_acp[indicatrice_pos,]-seuil_marg),
                                 FUN=fonction_trajectoire_positive,MARGIN = 1)==TRUE)
    indicatrice_pos<-indicatrice_pos[indicatrice_pos2]
    # #Truncate method
    # #If curve below the threshold, give the value of the threshold. 
    # Inds_diff<-Nouvelles_obs_P_acp<seuil_marg
    # Vect_transf<-c(Nouvelles_obs_P_acp)
    # print(summary(Vect_transf))
    # Vect_transf[Inds_diff]<-seuil_marg+10^(-10)
    # print(summary(Vect_transf))
    # Nouvelles_obs_P_acp<-matrix(Vect_transf,ncol=ncol(Nouvelles_obs_P_acp))
    # indicatrice_pos<-which(apply(X=(Nouvelles_obs_P_acp-seuil_marg),
    #                              FUN=fonction_trajectoire_positive,MARGIN = 1)==TRUE)
  }
  else{
    print("iciP")
    indicatrice_pos<-which(apply(X=Nouvelles_obs_P_acp-1,
                                 FUN=fonction_trajectoire_positive,MARGIN = 1)==TRUE)
    
  }
  Nouvelles_obs_P_acp<-Nouvelles_obs_P_acp[indicatrice_pos,]
  f_copy<-FORME_ACP[indicatrice_pos,]
  Ecart<-M-length(indicatrice_pos)
  forme_pr_analyse<-f_copy
  coord_pr_analyse<-Matrice_coords_ALL[indicatrice_pos,]
  while (Ecart!=0){
    if(NB_dim_ACP==2){
      Real_copules<-BiCopSim(N=Ecart,family =Famille_choisie$family,par = Famille_choisie$par,par2 =Famille_choisie$par2,check.pars = TRUE)
    }
    else{
      Real_copules<-VineCopula::RVineSim(N=Ecart,RVM = Famille_choisie)
    }
    Matrice_coords_ecart<-matrix(NA,ncol = NB_dim_ACP,nrow = Ecart)
    for(l in c(1:NB_dim_ACP)){
      if(Ecart==1){
        Matrice_coords_ecart[,l]<-quantile(Coordonnees[,l],probs=Real_copules[,l])
      }
      else{
        Matrice_coords_ecart[,l]<-quantile(Coordonnees[,l],probs=Real_copules[l])
      }
    }
    FORME_ecart<-function_reconstitution_trajectory_std(Vector_coords = Matrice_coords_ecart,
                                                        Base_functions_p  = F_propres,NB_dim = NB_dim_ACP,
                                                         mu_t =esp_temps,sd_t = sd_temps)
    #FORME_ecart<-exp(FORME_ecart)
    L2_ecart<-apply(X =FORME_ecart,MARGIN = 1,FUN = calcul_norme_L2)
    FORME_ecart<-t(t(FORME_ecart)%*%diag(L2_ecart^(-1)))
    
    valeurs_unif_e<-runif(Ecart)
    # Nouvelles echelles suivent une Pareto(1).
    reals_P_ecart<-sapply(valeurs_unif_e,Generation_Pareto_std)
    ECHELLE_ecart<-reals_P_ecart*v_seuil
    
    Nouvelles_obs_ecart<-t(t(FORME_ecart)%*%diag(ECHELLE_ecart))
    if(opt_Frech==TRUE){
      indicatrice_pos<-which(apply(X=Nouvelles_obs_ecart,
                                   FUN=fonction_trajectoire_positive,MARGIN = 1)==TRUE)
      Unif<-exp(-Nouvelles_obs_ecart[indicatrice_pos,]^(-1))
      PTO<-(1-Unif)^(-1)
      #seuil_marg<-1+10^(-5)
      seuil_marg<-log(10^(5))^(-1)
      indicatrice_pos2<-which(apply(X=(Nouvelles_obs_ecart[indicatrice_pos,]-seuil_marg),
                                    FUN=fonction_trajectoire_positive,MARGIN = 1)==TRUE)
    }
    else{
      indicatrice_pos<-which(apply(X=Nouvelles_obs_ecart-1,
                                   FUN=fonction_trajectoire_positive,MARGIN = 1)==TRUE)
    }
    if(length(indicatrice_pos)!=0){
      f_copy<-FORME_ecart[indicatrice_pos,]
      forme_pr_analyse<-rbind(forme_pr_analyse,f_copy)
      coord_pr_analyse<-rbind(coord_pr_analyse,
                              Matrice_coords_ecart[indicatrice_pos,])
      Nouvelles_obs_P_acp<-rbind(Nouvelles_obs_P_acp,
                                 Nouvelles_obs_ecart[indicatrice_pos,])
      Ecart<-Ecart-length(indicatrice_pos)
    }
    
  }
  L2_realisations<-apply(X = Nouvelles_obs_P_acp,MARGIN = 1,FUN = calcul_norme_L2)
  INDICES_ACP<-1:nrow(Nouvelles_obs_P_acp)
  Angle_realisations<-t(t(Nouvelles_obs_P_acp)%*%diag(L2_realisations^(-1)))
  write.csv(x = Nouvelles_obs_P_acp,
            file=paste0("residus_clust/Pareto/",nom_variable,"_simul_Pareto.csv"))

  # Comparaison_coordonnees -------------------------------------------------
  write.csv(x=vecteur_pareto,file=paste0("residus_clust/Pareto/",nom_variable,"_obs_ech_Pareto.csv"))
  #forme_pr_analyse_orig<-log(forme_pr_analyse+0.1)
  
  Min_y<-min(c(apply(forme_pr_analyse,MARGIN = 2,FUN = min),
               apply(Comp_angulaire,MARGIN = 2,FUN = min)))
  Max_y<-max(c(apply(forme_pr_analyse,MARGIN = 2,FUN = max),
               apply(Comp_angulaire,MARGIN = 2,FUN = max)))
  Time<-(c(1:37)-19)*(1/6)
  xlab_chosen<-"Time (hour) with respect to tidal peak"
  matplot(Time,t(Comp_angulaire),type="l",ylim=c(Min_y,Max_y),xlab=" ",
          ylab="value (-)",cex.lab=1.5)
  #mtext("Angles observed in the data")
  matplot(Time,t(forme_pr_analyse),type="l",ylim=c(Min_y,Max_y), 
          xlab=NA,ylab=NA,cex.lab=1.5)
  nom_graph<-ifelse(nom_variable=="Surcote","S",
                    nom_variable)
  #mtext(paste0("Simulated shape for ",nom_graph," with J=",NB_dim_ACP))
  mtext(xlab_chosen, side = 3, line = -22, outer = TRUE,
        cex=1.5)
  
  
  # travail norme -----------------------------------------------------------
  
  norme_theta_simul<-apply(X = forme_pr_analyse,MARGIN = 1,
                           FUN = calcul_norme_L2)
  par(mfrow=c(1,1))
  plot(density(norme_theta_simul),
       main=expression(paste0("Densité de la norme L2 de ",theta," simulé")))
  abline(v=1,col="red")
  write.csv(file=paste0("residus_clust/ACP_results/Theta_Vine_simul_",nom_variable,".csv"),
            x=coord_pr_analyse)
  write.csv(file=paste0("residus_clust/ACP_results/evolution_eigen_value_",nom_variable,".csv"),
            x=val_lambda)
  

  if(std_donnees==TRUE){
    Angle_realisations_std<-scale(Angle_realisations,center = esp_temps,
                                  scale = sd_temps)
    Conversion_realisations<-Angle_realisations_std%*%F_propres[,c(1:NB_dim_ACP)]
  }
  else{
    Conversion_realisations<-Angle_realisations%*%F_propres[,c(1:NB_dim_ACP)]
    
  }
  Matrice_test_nulles<-matrix(nrow=2,ncol=2*NB_dim_ACP)
  Z<-1
  for(j in c(1:NB_dim_ACP)){
    KS_j<-ks.test(x = Conversion_realisations[,j],
                  y=Coordonnees[,j])$p.value
    MW_j<-wilcox.test(x = Conversion_realisations[,j],
                      y=Coordonnees[,j])$p.value
    Matrice_test_nulles[,Z]<-c(KS_j,MW_j)
    Z<-Z+1
  }

  # Export coords pour echelle Pareto ---------------------------------------
  write.csv(file=paste0("residus_clust/ACP_results/",nom_variable,"_coords_obs.csv")
            ,x=Coordonnees)
  write.csv(file=paste0("residus_clust/ACP_results/",nom_variable,"_coords_simul.csv")
            ,x=Conversion_realisations)
  # Faire avec opt_Frech. Reconversion en Pareto pour 
  # la reconversion. 
  Avant_transf<-Nouvelles_obs_P_acp
  if(opt_Frech==TRUE){
    Unif<-exp((-1)*Nouvelles_obs_P_acp^(-1))
    Nouvelles_obs_P_acp<-(1-Unif)^(-1)
  }
  NO_ACP<-lapply(INDICES_ACP,FUN = fnct_select_colonne,df=Nouvelles_obs_P_acp)

  # Reconversion dans la bonne échelle --------------------------------------
  Variables_reconversion_ACP<-parLapply(cl = coeurs,NO_ACP,function_reconversion_Pareto,K=K,list_evt=L_EVT)
  Variables_reconversion_ACP<-do.call(rbind.data.frame,
                                      Variables_reconversion_ACP)
  
  colnames(Variables_reconversion_ACP)<-1:ncol(Variables_reconversion_ACP)
  # BResnick ----------------------------------------------------------------
  weigthFun <- function(x, u){
    return(x * (1 - exp(-((calcul_norme_L2(x/ u) - 1)))))

  }
  #Define partial derivative of weighting function
  dWeigthFun <- function(x, u){
    return((1 - exp(-((calcul_norme_L2(x/ u) - 1)))) +(x^2/u)*d_cste_norme_L2(x/u)*exp(-((calcul_norme_L2(x/ u) - 1))))

  }
  LOC<-expand.grid(1:37/37)
  vecteur_temps<-c(1:37)/37
  Excedents<-as.list(as.data.frame(t(Exts)))
  if(type_vario=="pow"){
    
    f_h_param<-function(h,parametre){
      return((abs(h)/parametre[1])^(parametre[2]))
    }
    objectiveFunction = function(parameter, exceedances, loc, vario, weigthFun, dWeigthFun, threshold){
      if(parameter[1]<0|parameter[2]<0|parameter[2]>=2){return(1e50)}
      varioModel <- function(h){
        vario(h = h,param = parameter)
      }
      mvPot::scoreEstimation(obs=exceedances, loc=loc, vario=varioModel, weigthFun, dWeigthFun, u = threshold)
    }
  }
  if(type_vario=="cubic"){
    f_h_param<-function(h,parametre){
      return(parametre[2]*(7*(abs(h)/parametre[1])^(2)-(35/4)*(abs(h)/parametre[1])^(3)
                           +(7/2)*(abs(h)/parametre[1])^(5)-(3/4)*(abs(h)/parametre[1])^(7)))
    }
    objectiveFunction = function(parameter, exceedances, loc, vario, weigthFun, dWeigthFun, threshold){
      if(parameter[1]<0|parameter[2]<0){return(1e50)}
      varioModel <- function(h){
        vario(h = h,param = parameter)
      }
      mvPot::scoreEstimation(obs=exceedances, loc=loc, vario=varioModel, weigthFun, dWeigthFun, u = threshold)
    }
  }
  if(type_vario=="exp"){
    f_h_param<-function(h,parametre){
      return(parametre[2]*(1-exp(-(abs(h)/parametre[1]))))
    }
    objectiveFunction = function(parameter, exceedances, loc, vario, weigthFun, dWeigthFun, threshold){
      if(parameter[2]<0|parameter[1]<0){return(1e50)}
      varioModel <- function(h){
        vario(h = h,param = parameter)
      }
      mvPot::scoreEstimation(obs=exceedances, loc=loc, vario=varioModel, weigthFun, dWeigthFun, u = threshold)
    }
  }
  if(type_vario=="powexp"){
    f_h_param<-function(h,parametre){
      return(parametre[3]*(1-exp(-(abs(h)/parametre[1])^(parametre[2]))))
    }
    objectiveFunction = function(parameter, exceedances, loc, vario, weigthFun, dWeigthFun, threshold){
      if(parameter[2]<0|parameter[1]<0|parameter[3]<0){return(1e50)}
      varioModel <- function(h){
        vario(h = h,param = parameter)
      }
      mvPot::scoreEstimation(obs=exceedances, loc=loc, vario=varioModel, weigthFun, dWeigthFun, u = threshold)
    }
  }
  if(type_vario=="gamma"){
    f_h_param<-function(h,parametre){
      return(parametre[3]*(1-(1+(abs(h)/parametre[1]))^(-parametre[2])))
    }
    objectiveFunction = function(parameter, exceedances, loc, vario, weigthFun, dWeigthFun, threshold){
      if(parameter[3]<0|parameter[2]<0|parameter[1]<0){return(1e50)}
      varioModel <- function(h){
        vario(h = h,param = parameter)
      }
      mvPot::scoreEstimation(obs=exceedances, loc=loc, vario=varioModel, weigthFun, dWeigthFun, u = threshold)
    }
  }
  
  if(type_vario=="gaussian"){
    
    f_h_param<-function(h,parametre){
      return(parametre[3]*(1-exp(-3*(abs(h)/parametre[1])^(2))))
    }
    objectiveFunction = function(parameter, exceedances, loc, vario, weigthFun, dWeigthFun, threshold){
      if(parameter[3]<0){return(1e50)}
      varioModel <- function(h){
        vario(h = h,param = parameter)
      }
      mvPot::scoreEstimation(obs=exceedances, loc=loc, vario=varioModel, weigthFun, dWeigthFun, u = threshold)
    }
  }
  if(type_vario=="cubic"|type_vario=="exp"|type_vario=="pow"){
    # Estimate the parameter by optimization of the objective function
    est <- optim(par =c(1,1),fn = objectiveFunction,
                 exceedances = Excedents,loc = LOC,
                 vario = f_h_param,weigthFun = weigthFun,
                 dWeigthFun = dWeigthFun,threshold =v_seuil,
                 control = list(maxit = 4000))
  }
  else{
    est <- optim(par =c(1,1,1),fn = objectiveFunction,
                 exceedances = Excedents,loc = LOC,
                 vario = f_h_param,weigthFun = weigthFun,
                 dWeigthFun = dWeigthFun,threshold =v_seuil,
                 control = list(maxit = 4000))
  }
  PARAMS_opt<-est$par
  print("params found")
  print(PARAMS_opt)
  f_variogram<-function(s,t){
    h<-abs(t-s)
    return(f_h_param(h,parametre = PARAMS_opt))
  }
  f_sigma<-function(h){
    return(f_h_param(h,parametre = PARAMS_opt))
  }
  # Extrémogramme empiriqiue ------------------------------------------------
  
  L<-ncol(Exts)
  z<-1

  # ANGLES<-t(Trajectoire_log_norm(Vecteur_temps=vecteur_temps,f_variogram=f_variogram,nb=10*M,alpha=1,f_sigma = f_sigma))
  # variogram_model_score<-function(h){
  #   return(f_h_param(h,parametre = PARAMS_opt))
  # }
  # # Simuler la composante de forme ------------------------------------------
  # FORME_BR<-Procedure_MHastings(echantillons_log_norm = ANGLES,Longueur_echantillon = M)
  # FORME_BR<-t(do.call(cbind.data.frame,FORME_BR))
  # colnames(FORME_BR)<-1:ncol(FORME_BR)
  # indicatrice_posBR<-which(apply(X=FORME_BR,FUN=fonction_trajectoire_positive,MARGIN = 1)==TRUE)
  # valeurs_unif<-runif(length(indicatrice_posBR))
  # reals_P_std<-sapply(valeurs_unif,Generation_Pareto_std)
  # ECHELLE_BR<-reals_P_std*v_seuil
  # 
  # BR_inds<-t(t(FORME_BR[indicatrice_posBR,])%*%diag(ECHELLE_BR))
  # if(opt_Frech==TRUE){
  #   Unif_BR<-exp(-1/BR_inds)
  #   BR_inds<-(1-Unif_BR)^(-1)
  # }
  # write.csv(x=BR_inds,file=paste0("residus_clust/Pareto/",nom_variable,"_simul_BR_Pareto.csv"))
  # Indices_BR<-c(1:nrow(BR_inds))
  # NO_BR<-lapply(Indices_BR,FUN = fnct_select_colonne,df=BR_inds)
  # 
  # # Reconversion dans la bonne échelle--Brown-Resnick --------------------------------------
  # Variables_reconversion_BR<-lapply(NO_BR,function_reconversion_Pareto,K=K,list_evt=L_EVT)
  # Variables_reconversion_BR<-do.call(rbind.data.frame,Variables_reconversion_BR)
  # colnames(Variables_reconversion_BR)<-1:ncol(Variables_reconversion_BR)
  # 
  # Coordonnées_echelle_initiale --------------------------------------------
  L2_X<-apply(X = Donnes[Inds_Exts,],FUN=calcul_norme_L2,MARGIN=1)
  Comp_ANG_X<-t(t(Donnes[Inds_Exts,])%*%diag(L2_X^(-1)))
  PCA_X<-PCA(X =Comp_ANG_X,scale.unit =std_donnees,ncp =NB_dim_ACP,graph = FALSE)
  Coords_donnees_X<-PCA_X$ind$coord[,c(1:NB_dim_ACP)]
  
  L2_realisations<-apply(X = Variables_reconversion_ACP,FUN=calcul_norme_L2,MARGIN=1)
  Comp_ANG_Xsimul<-t(t(Variables_reconversion_ACP)%*%diag(L2_realisations^(-1)))
  if(std_donnees==TRUE){
    sdX_temps<-apply(X = Comp_ANG_X,MARGIN = 2,FUN = sd)
    espX_temps<-colMeans(Comp_ANG_X)
    Comp_ANG_Xsimul<-scale(Comp_ANG_Xsimul,center = espX_temps,scale = sdX_temps)
  }
  Coords_simul_X<-Comp_ANG_Xsimul%*%(PCA_X$svd$V[,c(1:NB_dim_ACP)])
  for(z in c(1:NB_dim_ACP)){
    KS_j<-ks.test(x = Coords_simul_X[,z],y=Coords_donnees_X[,z])$p.value
    MW_j<-wilcox.test(x = Coords_simul_X[,z],y=Coords_donnees_X[,z])$p.value
    Matrice_test_nulles[,Z]<-c(KS_j,MW_j)
    Z<-Z+1
  }
  # Export des simulations_resid --------------------------------------------
  Std<-ifelse(std_donnees,"avec_std","sans_std")
  if(type_entree=="residuals"){
    #write.csv(file=paste0("residus_clust/Simulations",nom_variable,"_",Std,"_BR",".csv"),x = Variables_reconversion_BR)
    write.csv(file=paste0("residus_clust/Simulations",nom_variable,"_",Std,".csv"),x = Variables_reconversion_ACP)
  }
  if(type_entree=="clust"){
    print("iciclust")
    Fichier<-paste0("../simulations_/univ/",nom_variable,"_clust.csv")
    Fichier2<-paste0("../simulations_/univ/",nom_variable,"_clust_BR.csv")
    write.csv(file=Fichier,x = Variables_reconversion_ACP)
    write.csv(file=Fichier2,x = Variables_reconversion_BR)
  }
  # else{
  #   Fichier<-paste0("../simulations_/univ/",nom_variable,"_simul_brutes.csv")
  #   Fichier2<-paste0("../simulations_/univ/",nom_variable,"_simul_brutes_BR.csv")
  #   write.csv(file=Fichier,x = Variables_reconversion_ACP)
  #   write.csv(file=Fichier2,x = Variables_reconversion_BR)
  # }
  
  # Enregistrement de l'indice des individus extsm -------------
  Data_ext<-Donnes[Inds_Exts,]
  if(type_entree=="clust"){
    Inds_Exts<-Inds_Exts[Vrais_indices]
  }
  Inds_Exts_mat<-cbind(Inds_Exts,c(1:length(Inds_Exts)))
  write.csv(Inds_Exts_mat,file = paste0("residus_clust/inds_extremes_donnees_",nom_variable,".csv"))
  Transf_for_return<-as.data.frame(Angle_realisations)
  colnames(Transf_for_return)<-c(1:ncol(Transf_for_return))
  #"Simulations_BR"=Variables_reconversion_BR,
  return(list("Simulations"=Variables_reconversion_ACP,"obs_exts"=Data_ext,
              "Mat_HO"=Matrice_test_nulles,
              "Angle_Pareto"=Transf_for_return,
              "Angle_Pareto_realite"=Comp_angulaire))
  
}

# Prendre en compte la standardisation ou non.  ---------------------------
type_entree<-"residuals"
nom_variable<-"Surcote"
STD_donnees<-TRUE
P_U<-rep(0.10,37)
PL2<-0.05
NB_scores<-3
opt_Frech<-TRUE
#10000
N_simul<-2000

type_V<-"powexp"

if(type_entree=="residuals"){
  l_donnees<-paste0("residus_clust/",nom_variable,"_residus.csv")
  racine_ML<-c("resultat_ML/residus/")
}
if(type_entree=="donnees_brutes"){
  repertoire<-"../ss_tend/HIVER/"
  l_donnees<-paste0(repertoire,nom_variable,"_ss_tend.csv")
  racine_ML<-c("resultat_ML/donnees_brutes/")
}
if(type_entree=="clust"){
  l_donnees<-paste0("residus_clust/clust_",nom_variable,".csv")
  
}
Exts<-Simulation_ACP_TS_EXT(type_vario = type_V,type_entree = type_entree,
                            type_donnees="hiver",coeurs = coeurs,
                            NB_dim_ACP = NB_scores,lien_donnees =l_donnees,
                            nom_variable = nom_variable,p_U = P_U,
                            p_L2 = PL2,n.dens =2^(14),opt_Frech =opt_Frech,
                            M=N_simul,std_donnees=STD_donnees)


Nb_ech<-nrow(Exts$obs_exts)
Nb_simul<-nrow(Exts$Simulations)
Base_ind_s<-c(1:Nb_simul)
Base_ind_r<-c(1:Nb_ech)
Indices_train_r<-sample(Base_ind_r,size = round(Nb_ech*0.7),replace = FALSE)
Indices_test_r<-which(!(Base_ind_r%in% Indices_train_r))

Indices_s<-sample(Base_ind_s,size =Nb_ech,replace = FALSE)
nb_train<-round(0.70*Nb_ech)
Indices_train_s<-Indices_s[c(1:nb_train)]
debut<-nb_train+1
Indices_test_s<-Indices_s[c(debut:Nb_ech)]

Exts$Mat_HO
# Outils ML.  -------------------------------------------------------------
source("../fonctions/fonctions_perfs_ML.R")
# Matrice de confusion selon le noyau utilisé -----------------------------
nom_var1<-nom_variable
matplot(t(Exts$Simulations),type="l")
ech_train<-list("simul"=Exts$Simulations[Indices_train_s,],
                "real"=Exts$obs_exts[Indices_train_r,])

ech_test<-list("simul"=Exts$Simulations[Indices_test_s,],
               "real"=Exts$obs_exts[Indices_test_r,])
Liste<-fnct_ML_realite_simulations_all_kernels(ech_test = ech_test,
                                               ech_train = ech_train,
                                               variable=nom_var1)
Liste_noyaux<-names(Liste)[c(1:4)]


# Modèle logistique--Avec glm ----------------------------------------------------------------
# # -----------------------------------------------------------------------
Liste_fonction_lien<-c("logit", "probit", "cloglog","log","cauchit")
LISTE_resultat_glm<-fnct_logistigue_glm_mult_links(ech_test = ech_test,
                                                   ech_train = ech_train,
                                                   liste_fncts = Liste_fonction_lien,
                                                   variable=nom_var1)
if(type_entree=="donnees_brutes"){

  # Si travail directement avec les données brutes --------------------------
  # Conserver en mémoire performance SVM ------------------------------------
  write.csv(file = paste0(racine_ML,"glm_perfs_",nom_variable,".csv"),x=LISTE_resultat_glm$performances_globales)
  write.csv(file = paste0(racine_ML,"svm_perfs_",nom_variable,".csv"),x=Liste$performances_globales)
}

Titre_graph<-paste0("Courbe ROC basee sur les ", type_entree," pour ",nom_variable)

AUC<-Plot_ROC_curves_mult_links(obj_link_GLM = LISTE_resultat_glm,titre_graphique = Titre_graph,
                           liste_cles =c("logit","probit","cauchit"),
                           racine_repertoire=racine_ML,
                           variable=nom_var1)

# Angle_echelle_originale -------------------------------------------------
NL2_simul<-apply(X=Exts$Simulations,MARGIN = 1,
                 FUN = calcul_norme_L2)
NL2_ext<-apply(X=Exts$obs_exts,MARGIN = 1,
               FUN = calcul_norme_L2)
Theta_ech_orig<-t(t(Exts$obs_exts)%*%diag(NL2_ext^(-1)))
Theta_simul_orig<-t(t(Exts$Simulations)%*%diag(NL2_simul^(-1)))

nom_theta<-paste0("Angle ",nom_variable)

ech_train<-list("simul"=Theta_simul_orig[Indices_train_s,],
                "real"=Theta_ech_orig[Indices_train_r,])

ech_test<-list("simul"=Theta_simul_orig[Indices_test_s,],
                "real"=Theta_ech_orig[Indices_test_r,])

SVM<-fnct_ML_realite_simulations_all_kernels(ech_test = ech_test,
                                             ech_train = ech_train,
                                        variable=nom_theta)
SVM$radial$resultat_conf_test
GLM<-fnct_logistigue_glm_mult_links(ech_test = ech_test,
                                    ech_train = ech_train,
                                    liste_fncts = Liste_fonction_lien,
                                    variable=nom_theta)
GLM$logit$resultat_conf_test
nom_L2<-paste0("Norme de ",nom_variable)


# Travail norme L2 --------------------------------------------------------

ech_train<-list("simul"=NL2_simul[Indices_train_s],
                "real"=NL2_ext[Indices_train_r])

ech_test<-list("simul"=NL2_simul[Indices_test_s],
               "real"=NL2_ext[Indices_test_r])

SVM_l2<-fnct_ML_realite_simulations_all_kernels(ech_test = ech_test,
                                                ech_train = ech_train,
                                            variable=nom_L2)
GLM_l2<-fnct_logistigue_glm_mult_links(ech_test = ech_test,
                                    ech_train = ech_train,
                                    liste_fncts = Liste_fonction_lien,
                                    variable=nom_L2)
GLM_l2$logit$resultat_conf_test
GLM_l2$probit$resultat_conf_test
GLM_l2$log$resultat_conf_test

for(kernel_ut in Liste_noyaux){
  x_vect<-rbind(Liste[[kernel_ut]]$resultat_conf_test,
                SVM_l2[[kernel_ut]]$resultat_conf_test,
                SVM[[kernel_ut]]$resultat_conf_test)
  if(!is.null(x_vect)){
    write.csv(paste0(racine_ML,"svm/confusion_test_",kernel_ut,"_",nom_variable,".csv"),x=x_vect)
  }
}
Titre_graph_l2<-paste0("Courbe ROC basee sur les ", type_entree," pour la norme L2 de ",nom_variable)
AUC_l2<-Plot_ROC_curves_mult_links(obj_link_GLM = GLM_l2,titre_graphique = Titre_graph_l2,
                                  liste_cles =Liste_fonction_lien,
                                  racine_repertoire=racine_ML,
                                  variable=paste0("norme L2 de ",nom_variable))

Titre_graph<-paste0("Courbe ROC basee sur les ", type_entree," pour angle de ",nom_variable)
AUC_2<-Plot_ROC_curves_mult_links(obj_link_GLM = GLM,titre_graphique = Titre_graph,
                           liste_cles =c("logit","probit","cauchit"),
                           racine_repertoire=racine_ML,
                           variable=paste0("angle ",nom_variable))
for(lien in Liste_fonction_lien){
  x_vect<-rbind(LISTE_resultat_glm[[lien]]$resultat_conf_test,
                GLM_l2[[lien]]$resultat_conf_test,
                GLM[[lien]]$resultat_conf_test)
  if(!is.null(x_vect)){
    
    write.csv(paste0(racine_ML,"glm/confusion_test_",lien,"_",nom_variable,".csv"),x=x_vect)
  }
}

# Comparaison des normes L2 -----------------------------------------------

norme_<-c(NL2_simul,NL2_ext)
origine<-c(rep("simul",length(NL2_simul)),rep("obs",length(NL2_ext)))
df<-cbind.data.frame(norme_,origine)
colnames(df)<-c("Norme","Origine")

# Mettre dans aes width ou alpha rajoute une légende supplémentaire
# sur le graphique. 
# Avec l'effet présent quand même. 

ggplot(data=df,aes(y=Norme,x=Origine,fill=Origine))+
  geom_violin(alpha=0.6)+
  geom_boxplot(width=0.3)+
  ylab("valeur")+
  ggtitle(paste0("Distribution de la norme L2 de ",nom_variable," selon l'origine (Residus)"))+
  labs(fill="Légende")

# Loi suivie par la norme L2 des simulations (vs) norme L2 donnees ----------------------------------------------
ACF_graph<-acf(NL2_simul,lag=500,plot = FALSE)
PACF_graph<-pacf(NL2_simul,lag=500,plot = FALSE)
plot(ACF_graph,ylim=c(-1,1),main="ACF série des normes simulées")
plot(PACF_graph,ylim=c(-1,1),main="PACF série des normes simulées")

# Exporter l'AUC --------------------------------------------------------------
print(AUC$grap)
print(AUC_2$grap)
print(AUC_l2$grap)

AUC_all<-rbind(AUC$AUC,AUC_2$AUC,AUC_l2$AUC)
AUC_all$methode<-type_entree
write.csv(x=AUC_all,file =paste0(racine_ML,"AUC_",type_entree,".csv"),row.names = FALSE)


# Fin code ----------------------------------------------------------------
# ------------------------------------------------------------------------
stopCluster(coeurs)


