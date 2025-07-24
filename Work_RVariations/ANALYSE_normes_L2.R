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

fonction_queue_lourde_L2<-function(coeurs,nom_variable,opt_tend,type_donnees,Min_u,Max_u){
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
  
  # ACF resultat ------------------------------------------------------------
  
  colnames(obs)<-c(1:length(colnames(obs)))
  VECTEUR_normeL2<-parApply(cl=coeurs,X = obs,MARGIN = 1,FUN = calcul_norme_L2)
  nb_obs<-750
  vect_y<-20:nb_obs
  k1<-sum(as.numeric(VECTEUR_normeL2>Min_u))
  k2<-sum(as.numeric(VECTEUR_normeL2>Max_u))
  ML_extRemes<-sapply(X =vect_y ,FUN = fonction_ML_extRemes,donnees=VECTEUR_normeL2)
  #-0.75 et 0.75 avant
  print(fonction_MLplot_resume(resultatML = ML_extRemes,vecteur_k = vect_y,nom_variable = paste0("la norme L2 de ",nom_variable),lims_Y =c(-0.75,0.75) ))

  # Hill --------------------------------------------------------------------
  boxplot(VECTEUR_normeL2)
  Hill_gamma<-sapply(X = vect_y,FUN = fonction_estimateur_hill,series_=VECTEUR_normeL2)
  plot(vect_y,Hill_gamma,type="l")
  
  # Moment ------------------------------------------------------------------
  Moments_gamma<-sapply(X = vect_y,FUN=fonction_estimateur_moment,series_=VECTEUR_normeL2)
  TITRE<-paste0("Estimateur du parametre de forme (moment) pour ",nom_variable," pour les donnees ", type_donnees)
  plot(vect_y,Moments_gamma,type="l",main=TITRE,ylab="Estimateur methode moments",xlab="Nombre d'excédents")
  abline(h=0,col="red")
  abline(v=k1,col="green")
  abline(v=k2,col="green")
  ML_gamma<-unlist(ML_extRemes[3,])
  data_Hill<-cbind.data.frame(vect_y,Hill_gamma)
  colnames(data_Hill)<-c("nombre_excedents","estimateur_Hill")
  GGHill<-ggplot(data = data_Hill,aes(x=nombre_excedents,y=estimateur_Hill))+
    geom_line()+
    ggtitle(paste0("Estimateur de Hill du paramètre de forme pour la norme L2 de Hs "))+
    labs(caption=paste0("k compris entre ",min(vect_y)," et ",max(vect_y)))
  print(GGHill)
  data_estimateur<-cbind(Moments_gamma,ML_gamma)
  Gam<-as.data.frame(melt(t(data_estimateur)))
  
  Gam$Var2<-vect_y[Gam$Var2]
  colnames(Gam)<-c("source","nombre_excedents","gamma_estime")
  GGothers<-ggplot(data=Gam,aes(x=nombre_excedents,y =gamma_estime,col=source,group=interaction(source)))+
    geom_line()+
    labs(colour="Legende")+
    ggtitle(paste0("Estimateur du paramètre de forme pour la norme L2 de ",nom_variable,""))+
    labs(caption=paste0("k compris entre ",min(vect_y)," et ",max(vect_y)))
  print(GGothers)
  
  # Choix du seuil ----------------------------------------------------------
  df_time_variable<-cbind.data.frame(1:length(VECTEUR_normeL2),VECTEUR_normeL2)
  colnames(df_time_variable)<-c("time","obs")
  
  # Excedents_dummy ---------------------------------------------------------
  ################
  seuil<-sort(VECTEUR_normeL2,decreasing = TRUE)[210]
  ML_dummy<-extRemes::fevd(x = VECTEUR_normeL2,threshold = seuil,type="GP")
  Params_dummy<-ML_dummy$results$par
  Gamma_ML<-as.numeric(Params_dummy[[2]])
  scale<-as.numeric(Params_dummy[[1]])
  excedents_dummy<-subset(VECTEUR_normeL2,VECTEUR_normeL2>seuil)
  print(Gamma_ML)
  QGPDpareto<-function(x){
    # si gamma est proche de 0. 
    if(abs(Gamma_ML)<10^{-4}){
      part1<-(-1)*log(1-x)*scale
      return(part1+seuil)
    }
    numerateur<-(1-x)^(-Gamma_ML)-1
    denom<-Gamma_ML
    y<-((numerateur/denom)*scale)+seuil
    return(y)
  }
  N_excedents<-length(excedents_dummy)

  # Première méthode (dummy) ------------------------------------------------
  simulationsGPD<-QGPDpareto(ppoints(n = N_excedents))
  extRemes::qqplot(excedents_dummy,simulationsGPD,xlab = "Quantiles empiriques",ylab="Quantiles théoriques")
  mtext(paste0("Comparaison de la norme L2 de ",nom_variable," en ",tolower(type_donnees)," avec une GPD(",round(seuil,1),",",round(scale,1),",",round(Gamma_ML,1),")"))
  KS_L2_dummy<-ks.test(excedents_dummy,extRemes::"pevd",scale=scale,shape=Gamma_ML,threshold=seuil,type="GP")
  AD_L2_dummy<-goftest::ad.test(x = excedents_dummy,null = extRemes::"pevd",
                                scale=scale,shape=Gamma_ML,threshold=seuil,
                                type="GP")
  print(KS_L2_dummy)
  print(AD_L2_dummy)
  # Travail fonction quantile -----------------------------------------------
  ############

  # (1) travail avec les quantiles ------------------------------------------
  
  p_normale<-sort(seq.int(from = 0.025,to = 0.20,length.out = 30)^(-1))
  fonction_GPD_exp_quantile<-function(scale,shape,seuil,p){
    part1<-p^(-shape)-1
    return((part1/shape)*scale+seuil)
  }
  Ninter<-1000
  simulations_for_int<-QGPDpareto(ppoints(n =Ninter))
  print(c(seuil,scale))
  valeurs_associees_theorique<-sapply(p_normale^(-1),FUN = fonction_GPD_exp_quantile,scale=scale,shape=Gamma_ML,seuil=seuil)
  valeurs_associees_empirique<-as.numeric(quantile(excedents_dummy,1-p_normale^(-1)))
  Quantiles_full_theor<-extRemes::qqplot(simulations_for_int,simulations_for_int,make.plot=FALSE)
  upper<-Quantiles_full_theor$qdata$upper
  lower<-Quantiles_full_theor$qdata$lower
  
  indices<-which((is.na(lower)==FALSE)&(is.na(upper)==FALSE))
  B_supKs<-as.numeric(quantile(upper[indices],1-p_normale^(-1)))
  B_infKs<-as.numeric(quantile(lower[indices],1-p_normale^(-1)))
  aspect_extreme<-cbind.data.frame(valeurs_associees_empirique,valeurs_associees_theorique)
  colnames(aspect_extreme)<-c("empirique","GPD")
  Melteur<-melt(t(aspect_extreme))
  print(quantile(VECTEUR_normeL2,0.85))
  print(quantile(VECTEUR_normeL2,0.98))
  colnames(Melteur)<-c("source","probabilite_depassement","valeur")
  Melteur$sup_KS<-B_supKs[Melteur$probabilite_depassement]
  Melteur$inf_KS<-B_infKs[Melteur$probabilite_depassement]
  Melteur$probabilite_depassement<-p_normale[Melteur$probabilite_depassement]
  Cols_<-c("empirique"="blue","GPD"="green","Percentile_5_95_KS"="darkblue")
  
  GGreturn_log<-ggplot(data=Melteur,aes(x =probabilite_depassement,y=valeur,color=source,group=interaction(source)))+
    geom_point()+
    geom_line()+
    scale_x_continuous(trans = "log10")+
    labs(caption=paste0("simulations_GPD_KS=",Ninter))+
    scale_color_manual(values=Cols_)+
    xlab("1/p")+
    geom_ribbon(data =Melteur,mapping = aes(ymin=inf_KS,ymax=sup_KS,col="Percentile_5_95_KS"),alpha=0.15,fill="grey", linetype = "dashed")+
    ggtitle(paste0("Comportement de la norme L2 de ",nom_variable," (repère log(x)) pour les excédents"))
  print(GGreturn_log)
  
  # GPD (queue lourde) ---------------------------------------------------
  valeurs_associees_Queue_lourde<-sapply(p_normale^(-1),FUN = fonction_GPD_exp_quantile,scale=scale,shape=1,seuil=seuil)
  plot(p_normale,valeurs_associees_Queue_lourde,log="x")
    
  # (2) Travail avec répartition --------------------------------------------
  vecteur<-c(excedents_dummy,simulations_for_int)
  Type<-c(rep("données",length(excedents_dummy)),rep("GPD",length(simulations_for_int)))
  df_for_densite<-cbind.data.frame(vecteur,Type)
  colnames(df_for_densite)<-c("valeur","source")
  Cols_<-c("données"="blue","GPD"="green")
  GG_dens<-ggplot(data=df_for_densite,aes(valeur,col=source,fill = source))+
    geom_density(alpha=0.1)+
    scale_color_manual(values = Cols_)+
    ggtitle(paste0("Comportement de la densité de la norme L2 de ",nom_variable," pour les excédents"))
  print(GG_dens)
  # Diagnostic du package pour affiner --------------------------------------
  ###########
  
  # travail sur l'ensemble de la norme --------------------------------------
  df_ensemble<-cbind.data.frame(VECTEUR_normeL2,rep("données",length(VECTEUR_normeL2)))
  colnames(df_ensemble)<-c("Norme_L2","origine")
  gillustr<-ggplot(data = df_ensemble,aes(x=Norme_L2,fill=origine))+
    geom_density()+
    ggtitle(paste0("Densité de la norme L2 pour ",nom_variable))+
    ylab("Densité")+
    scale_fill_manual(values = c("données"="blue"))+
    xlab("Valeur de la norme")+
    geom_vline(aes(xintercept = seuil,col="seuil"))+
    labs(colour="légende")
  print(gillustr)

  par(mfrow=c(2,2))
  valeurs<-POT::tcplot(VECTEUR_normeL2,ask = FALSE,u.range = c(Min_u,Max_u))
  title(paste("evolution scale/shape, variable=",nom_variable,"sans tendance",opt_tend,"pour",type_donnees))
  POT::mrlplot(data = VECTEUR_normeL2,u.range = c(Min_u,Max_u))
  # # With time.
  POT::diplot(data = df_time_variable,u.range = c(Min_u,Max_u),nt=200)
  q1<-mean(as.numeric(VECTEUR_normeL2<Min_u))
  q2<-mean(as.numeric(VECTEUR_normeL2<Max_u))
  
  qu_threshr<- quantile(VECTEUR_normeL2, probs = seq(q1, q2,length.out=10))
  Resultat_threshr<-threshr::ithresh(data=VECTEUR_normeL2,u_vec =qu_threshr)
  seuil_L2<-summary(Resultat_threshr)[3]
  Qopt<-summary(Resultat_threshr)[4]
  Prop<-(Qopt/100)
  Indice_opt_gamma<-round(length(VECTEUR_normeL2)*(1-Prop))
  print(paste0("Le nombre d'excédents est ",Indice_opt_gamma))
  par(mfrow=c(1,1))
  plot(vect_y,Moments_gamma,type="l",main=TITRE,ylab="Estimateur methode moments",xlab="Nombre d'excédents")
  abline(v=Indice_opt_gamma,col="green")
  
  # Utiliser les params de la méthode ---------------------------------------
  ###############
  seuil<-quantile(VECTEUR_normeL2,Prop)
  modele_evd<-extRemes::fevd(VECTEUR_normeL2,threshold = seuil, type = "GP")
  plot(modele_evd)
  Params<-modele_evd$results$par
  Gamma_ML<-as.numeric(Params[[2]])
  scale<-as.numeric(Params[[1]])
  excedents<-subset(VECTEUR_normeL2,VECTEUR_normeL2>seuil)

  # Comparaison graphique ---------------------------------------------------
  ##################
  
  N_excedents<-length(excedents)
  extRemes::qqplot(excedents,QGPDpareto(ppoints(n = N_excedents)),ylab = "Quantiles théoriques",xlab="Quantiles empiriques")
  mtext(paste0("Comparaison de la norme L2 de ",nom_variable," en ",tolower(type_donnees)," avec une GPD(",round(seuil,1),",",round(scale,1),",",round(Gamma_ML,1),")"))
  
  # Test statistique --------------------------------------------------------
  ##############
  KS_L2_X<-ks.test(excedents,extRemes::"pevd",scale=scale,shape=Gamma_ML,threshold=seuil,type="GP")
  AD_L2_X<-goftest::ad.test(x = excedents,null = extRemes::"pevd",scale=scale,shape=Gamma_ML,threshold=seuil,type="GP")

  return(list("KS"=c(KS_L2_dummy,KS_L2_X),"AD"=c(AD_L2_dummy,AD_L2_X)))
}

R<-fonction_queue_lourde_L2(coeurs = coeurs,nom_variable = "Hs",opt_tend = TRUE,type_donnees = "HIVER",Min_u=3.5,Max_u=5.4)
#fonction_queue_lourde_L2(coeurs = coeurs,nom_variable = "Surcote",opt_tend = TRUE,type_donnees = "HIVER",Min_u=12,Max_u=15)
#Surcote : 0.197, 0.29
#resultat_surcote<-fonction_queue_lourde_L2(coeurs = coeurs,nom_variable = "Surcote",opt_tend = TRUE,type_donnees = "HIVER",Min_u=0.197,Max_u=0.29)
#resultat_surcote

# Fermeture parallel ------------------------------------------------------
stopCluster(coeurs)