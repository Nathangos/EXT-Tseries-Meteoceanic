
# Charger les seuils obtenus avec la méthode threshr ----------------------
rm(list=ls())
par(mfrow=c(1,1))
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
source("fonctions_Pareto.R")
F_names<-ls()
nb_coeurs<-detectCores()-4
coeurs<-makeCluster(nb_coeurs)
clusterExport(coeurs,varlist=F_names)

fonction_queue_lourde_L2_TRANSF<-function(coeurs,nom,n.dens,opt_tend,Tau,type_donnees,simple_P,Min_u,Max_u,lims_y_gamma){
  
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
  print(nrow(obs))
  ## Table pour la tendance --------------------------------------------------
  nom_rech2<-paste0(repertoire,"RL_",nom,".csv")
  MATRICE_tend<-read.csv(nom_rech2)[,2:38]
  if (opt_tend==FALSE){
    obs<-obs+MATRICE_tend
  }
  p_u<-c()
  for(j in c(1:ncol(obs))){
    prop<-mean(as.numeric(obs[,j]>Tau[j]))
    p_u<-c(p_u,prop)
  }
  ## Approcher le poids mis dans la queue de distribution.

  print("resume Donnees")
  print(sum(colMeans(obs)))
  colnames(obs)<-c(1:length(colnames(obs)))
  ## Preciser si la tendance est utilisee.
  Tendance_ounon<-ifelse(opt_tend,"On enlève la tendance","On conserve la tendance")
  print(Tendance_ounon)
  boxplot(obs,main=paste("Boxplot de la variable",nom))
  Vecteur_DIM<-1:ncol(obs)
  ## Donner des marginales uniformes sur (-1,0).
  Resultat_p<-as.data.frame(t(sapply(Vecteur_DIM,f_marginales_all_Pareto,donnees=obs,p_u=p_u,n.dens=n.dens,donnees_a_transformer=obs)))
  K_MATRICE_DENS<-Resultat_p$kernel_dens
  Pareto<-Resultat_p$obs
  Pareto<-do.call(cbind.data.frame,Pareto)
  colnames(Pareto)<-c(1:ncol(Pareto))
  P_valeur<-c()
  P_valeur_AD<-c()
  print("Poids mis")
  print(p_u)
  P_valeur_AD_excedent_GPD<-c()
  P_valeur_KS_excedent_GPD<-c()
  vect_Gamma<-c()
  
  # Observer le gamma en chaque temps ---------------------------------------
  for (t in 1:ncol(Pareto)){
    variable_t<-Pareto[,t]
    TEST_ks_pareto<-ks.test(x = variable_t,extRemes::"pevd",threshold=1,scale=1,shape=1,type="GP")
    P_valeur<-c(P_valeur,TEST_ks_pareto$p.value)
    TEST_AD_pareto<-goftest::ad.test(x = variable_t,extRemes::"pevd",threshold=1,scale=1,shape=1,type="GP")
    P_valeur_AD<-c(P_valeur_AD,TEST_AD_pareto$p.value)
    vecty<-c(50:(length(variable_t)/2))
    
    # Moment ------------------------------------------------------------------
    p_ut<-p_u[t]
    # ML_loi excedents --------------------------------------------------------
    variable_ech_original<-obs[,t]
    seuil_t<-quantile(variable_ech_original,1-p_ut)
    kernel_dens<-density(x=variable_ech_original,n=n.dens)
    p_u_t<-as.numeric(mean(variable_ech_original>seuil_t))
    indice<-which.max(seuil_t-kernel_dens$x<0)-1
    valeur_fonc_u<-as.numeric(kernel_dens$y[indice])
    sigma_t<-(p_u_t/valeur_fonc_u)
    #indice pour k. 
    gamma_t<-get.tail.index.estimate(variable_ech_original,1-p_u_t,sigma_t)
    queue_distrib<-subset(variable_ech_original,variable_ech_original>seuil_t)

    QGPDpareto<-function(x){
      # si gamma est proche de 0. 
      if(abs(gamma_t)<10^{-4}){
        part1<-(-1)*log(1-x)*sigma_t
        return(part1+seuil_t)
      }
      numerateur<-(1-x)^(-gamma_t)-1
      denom<-gamma_t
      y<-((numerateur/denom)*sigma_t)+seuil_t
      return(y)
    }
    N_excedents<-length(queue_distrib)
    qqplot(x=queue_distrib,y=QGPDpareto(ppoints(n = N_excedents)),ylab = "Quantiles théoriques",xlab="Quantiles empiriques")
    # qqline(y=queue_distrib,datax=FALSE,distribution=QGPDpareto)
    mtext(paste0("Comparaison avec une GPD (excédents) t=",t,", variable=",nom," en ",tolower(type_donnees)))
    
    
    M<-length(variable_t)
    qqplot(x=variable_t,y=Qpareto(ppoints(n = M)),ylab = "Quantiles théoriques",xlab="Quantiles empiriques")
    mtext(paste0("Comparaison avec une Pareto t=",t,", variable=T(",nom,") en ",tolower(type_donnees)))
    vect_Gamma<-c(vect_Gamma,gamma_t)
    TEST_AD_pareto_excedents<-goftest::ad.test(queue_distrib,extRemes::"pevd",threshold=seuil_t,scale=sigma_t,shape=gamma_t,type="GP")
    TEST_KS_pareto_excedents<-ks.test(queue_distrib,extRemes::"pevd",threshold=seuil_t,scale=sigma_t,shape=gamma_t,type="GP")
    P_valeur_AD_excedent_GPD<-c(P_valeur_AD_excedent_GPD,TEST_AD_pareto_excedents$p.value)
    P_valeur_KS_excedent_GPD<-c(P_valeur_KS_excedent_GPD,TEST_KS_pareto_excedents$p.value)
  }
  vect_temps<-c(1:37)
  saison<-tolower(type_donnees)
  print(vect_Gamma)
  plot(vect_temps,vect_Gamma,main=paste0("\u03B3", "(t) estime pour ",nom," en ",saison),xlab="t",ylab="Estimateur",type="o")
  df_test_excedent_GPD<-cbind.data.frame(P_valeur_KS_excedent_GPD,P_valeur_AD_excedent_GPD)
  colnames(df_test_excedent_GPD)<-c("KS","AD")
  GG_TEST_GPD<-ggplot2::ggplot(data = df_test_excedent_GPD,aes(x=1:37,y=KS,col="KS"))+
    geom_point()+
    geom_line(aes(y=KS,col="KS"))+
    geom_point(aes(y=AD,col="AD"))+
    geom_line(aes(y=AD,col="AD"))+
    ylim(c(0,1))+
    ggtitle(paste0("p valeur des tests ",nom,"(t)|",nom,"(t)>u(t)~GPD(u(t)",",","\u03C3","(t),","\u03B3","(t))"))+
    xlab(label = "t")+
    ylab(label="p valeur")+
    geom_hline(yintercept = 0.05,show.legend = TRUE)+
    scale_color_manual("Test",values = c("red","blue"))+
    labs(caption = paste0("Valeur du seuil moyen arrondie =",round(mean(Tau),2),", n=",length(variable_t),", n.dens=",n.dens))
  print(GG_TEST_GPD)
  
  df_test_P<-cbind.data.frame(P_valeur,P_valeur_AD)
  colnames(df_test_P)<-c("KS","AD")
  GG_TEST<-ggplot2::ggplot(data = df_test_P,aes(x=1:37,y=KS,col="KS"))+
    geom_point()+
    geom_line(aes(y=KS,col="KS"))+
    geom_point(aes(y=AD,col="AD"))+
    geom_line(aes(y=AD,col="AD"))+
    ylim(c(0,1))+
    ggtitle(paste0("p valeur des tests T(",nom,")(t)~Pareto(1)"))+
    xlab(label = "t")+
    ylab(label="p valeur")+
    geom_hline(yintercept = 0.05,show.legend = TRUE)+
    scale_color_manual("Test",values = c("red","blue"))+
    labs(caption = paste0("Valeur du seuil moyen arrondie =",round(mean(Tau),2),", n=",length(variable_t),", n.dens=",n.dens))
  print(GG_TEST)
  # Observer l'évaluation de la p_valeur ------------------------------------
  
  ## Calcul du vecteur des normes inverses.
  vecteur_Normes_Transf_Pareto<-apply(X=Pareto,MARGIN = 1,FUN = calcul_norme_L2)
  vecteur_somme<-apply(X = Pareto,MARGIN = 1,FUN = sum)
  nb_obs<-length(vecteur_Normes_Transf_Pareto)/3
  print("quantile pour le maximum")
  vect_y<-20:nb_obs
  ML_extRemes<-parSapply(cl = coeurs,X =vect_y ,FUN = fonction_ML_extRemes,donnees=vecteur_Normes_Transf_Pareto)
  print(fonction_MLplot_resume(resultatML = ML_extRemes,vecteur_k = vect_y,nom_variable = paste0("la norme L2 de T(",nom,")"),lims_Y=lims_y_gamma))
  # Pour la somme -----------------------------------------------------------
  ML_extRemes_somme<-parSapply(cl = coeurs,X =vect_y ,FUN = fonction_ML_extRemes,donnees=vecteur_somme)
  print(fonction_MLplot_resume(resultatML = ML_extRemes_somme,vecteur_k = vect_y,nom_variable = paste0("la norme L2 de T(",nom,")"),lims_Y=lims_y_gamma))
  
  # Hill --------------------------------------------------------------------
  Hill_gamma<-sapply(X = vect_y,FUN = fonction_estimateur_hill,series_=vecteur_Normes_Transf_Pareto)
  TITRE<-paste0("Estimateur du parametre de forme (Hill) pour T(",nom,") pour les donnees ", type_donnees)
  plot(vect_y,Hill_gamma,type="l",main=TITRE,ylab="Estimateur methode Hill",xlab="Nombre d'excédents")
  
  # Moment ------------------------------------------------------------------
  k1<-sum(as.numeric(vecteur_Normes_Transf_Pareto>Min_u))
  k2<-sum(as.numeric(vecteur_Normes_Transf_Pareto>Max_u))
  Moments_gamma<-sapply(X = vect_y,FUN=fonction_estimateur_moment,series_=vecteur_Normes_Transf_Pareto)
  TITRE<-paste0("Estimateur du parametre de forme (moment) pour T(",nom,") pour les donnees ", type_donnees)
  plot(vect_y,Moments_gamma,type="l",main=TITRE,ylab="Estimateur methode moments",xlab="Nombre d'excédents")
  abline(h=1,col="red")
  abline(v=k1,col="green")
  abline(v=k2,col="green")
  ML_gamma<-unlist(ML_extRemes[3,])
  data_estimateur<-cbind(Moments_gamma,Hill_gamma,ML_gamma)
  Gam<-as.data.frame(melt(t(data_estimateur)))

  Gam$Var2<-vect_y[Gam$Var2]
  colnames(Gam)<-c("source","nombre_excedents","gamma_estime")
  GGothers<-ggplot(data=Gam,aes(x=nombre_excedents,y =gamma_estime,col=source,group=interaction(source)))+
    geom_line()+
    labs(colour="Legende",caption=paste0("k compris entre ",min(vect_y)," et ",max(vect_y)))+
    ggtitle(paste0("Estimateur du paramètre de forme pour la norme L2 de T(",nom,")"))
  
  print(GGothers)
  # Choix du seuil ----------------------------------------------------------
  df_time_variable<-cbind.data.frame(1:length(vecteur_Normes_Transf_Pareto),vecteur_Normes_Transf_Pareto)
  colnames(df_time_variable)<-c("time","obs")
  par(mfrow=c(2,2))
  Sous_pop<-subset(vecteur_Normes_Transf_Pareto,vecteur_Normes_Transf_Pareto<max(vecteur_Normes_Transf_Pareto))
  plot(density(Sous_pop))
  U1<-quantile(vecteur_Normes_Transf_Pareto,Min_u)
  U2<-quantile(vecteur_Normes_Transf_Pareto,Max_u)
  valeurs<-POT::tcplot(vecteur_Normes_Transf_Pareto,ask = FALSE,u.range = c(U1,U2))
  title(paste0("Analyse de la norme L2 de T(",nom,")"),line = -3, outer = TRUE)
  POT::mrlplot(data = vecteur_Normes_Transf_Pareto,u.range = c(U1,U2))
  # # With time.
  Diplot_threshr<-POT::diplot(data = df_time_variable,u.range = c(U1,U2),nt=50)
  qu_threshr<- quantile(vecteur_Normes_Transf_Pareto, probs = seq(Min_u, Max_u,length.out=50))
  Resultat_threshr<-threshr::ithresh(vecteur_Normes_Transf_Pareto,u_vec =qu_threshr)
  par(mfrow=c(1,1))
  plot(Resultat_threshr)
  seuil_NORMEL2<-summary(Resultat_threshr)[3]
  Prop<-mean(as.numeric(vecteur_Normes_Transf_Pareto<seuil_NORMEL2))
  # Test Hypothese pareto pour excedents ------------------------------------
  
  vecteur_excedents<-vecteur_Normes_Transf_Pareto[which(vecteur_Normes_Transf_Pareto>seuil_NORMEL2)]
  M<-length(vecteur_excedents)
  qqplot(x=vecteur_excedents/seuil_NORMEL2,y=Qpareto(ppoints(n = M)),xlab = "Quantiles empiriques",ylab="Quantiles théoriques")
  mtext(paste0("Qqplot Pareto pour la (norme L2/seuil) de ",nom," en ",tolower(type_donnees)," (excédents)"))
  print(Prop)
  
  plot(density(vecteur_excedents),main="Densité norme L2 valeurs exts")
  TEST_ks_pareto<-ks.test(x = vecteur_excedents,extRemes::"pevd",threshold=seuil_NORMEL2,scale=seuil_NORMEL2,shape=1,type="GP")
  P_valeur<-TEST_ks_pareto$p.value
  TEST_AD_pareto<-goftest::ad.test(x = vecteur_excedents,extRemes::"pevd",
                                   threshold=seuil_NORMEL2,scale=seuil_NORMEL2,shape=1,
                                   type="GP")
  P_valeur_AD<-TEST_AD_pareto$p.value
  print(c(P_valeur,P_valeur_AD))
} 


fonct_lancement<-function(nom_variable,type_donnees,Min_u,Max_u,N,type_seuils){
  e_<-read.table(paste0("SEUILS_THRESHR_",nom_variable,"_",type_donnees,"_",type_seuils,".txt"))
  SEUILS_THRESHR<-e_$x
  fonction_queue_lourde_L2_TRANSF(coeurs = coeurs,nom = nom_variable,n.dens = N,opt_tend = TRUE,Tau =cummean(SEUILS_THRESHR),type_donnees = type_donnees,Max_u=Max_u,Min_u=Min_u,simple_P = NA,lims_y_gamma=c(0,2.3))
}
resultat<-fonct_lancement("Hs","HIVER",Min_u=0.85,Max_u=0.98,N=2^(13),type_seuils="POT")

# Fin parallel -----------------------------------------------------------
stopCluster(coeurs)