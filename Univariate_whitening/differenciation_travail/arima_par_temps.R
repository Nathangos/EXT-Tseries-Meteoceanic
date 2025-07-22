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
require(tseries)

### Source des fonctions. ####
source("fonctions/fonctions.R")
F_names<-ls()
nb_coeurs<-detectCores()-4
coeurs<-makeCluster(nb_coeurs)
clusterExport(coeurs,varlist=F_names)

dates_Johanna<-c("08/03/2008","09/03/2008","10/03/2008")
date_indu_ext_plus<-"17/12/1989"


fonction_arima_variable<-function(coeurs,nom_variable,opt_tend,type_donnees,p,q,P,D,Q,S,pas_clust,C_thchosen,C_Falk){
  repertoire<-"ss_tend/"
  if(type_donnees=="HIVER"){
    repertoire<-paste0(repertoire,"HIVER/")
  }
  if(type_donnees=="HORS_HIVER"){
    repertoire<-paste0(repertoire,"HORS_HIVER/")
  }
  nom_recherche<-paste0(repertoire,nom_variable,"_ss_tend.csv")
  obs<-read.csv(nom_recherche)[,2:38]
  if(nom_variable=="Surcote"){
    nom_variable<-"S"
  }
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
  # Nombre individus --------------------------------------------------------
  print(nrow(obs))
  
  # ACF resultat ------------------------------------------------------------
  B<-melt(t(obs))
  colnames(B)<-c("Temps","Individu","Mesure")
  B$Temps<-as.factor(B$Temps)
  TITRE<-paste0("Densité avant diff en chaque temps pour ",nom_variable, " en hiver")
  PV<-ggplot2::ggplot(data = B,aes(Temps,Mesure,fill=Temps))+
    geom_violin()+
    geom_boxplot(width=0.2)+
    guides(fill = "none")+
    ggtitle(TITRE)+
    ylab("Densité")
  
  colnames(obs)<-c(1:length(colnames(obs)))
  VECTEUR_normeL2<-parApply(cl=coeurs,X = obs,MARGIN = 1,FUN = calcul_norme_L2)
  

  # Choix seuil avant transformation ----------------------------------------
  Q1_brut<-quantile(VECTEUR_normeL2,0.1)
  Q2_brut<-quantile(VECTEUR_normeL2,0.95)
  qu_threshr_brut<- quantile(VECTEUR_normeL2, probs = seq(0.1, 0.95,length.out=10))
  #Resultat_seuil_brut<-threshr::ithresh(data =VECTEUR_normeL2,u_vec =qu_threshr_brut)
  
  par(mfrow=c(1,2))
  df_time_variable<-cbind.data.frame(1:length(VECTEUR_normeL2), VECTEUR_normeL2)
  colnames(df_time_variable)<-c("time","obs")
  dates_prises<-read.csv(file="ss_tend/HIVER/dates_prises.csv")[,2]
  # Convertir en nombre de jours-3 jours =une tempête.
  df_time_variable$time<-dates_prises
  df_time_variable$time<-lubridate::decimal_date(as.POSIXct(df_time_variable$time, format="%d/%m/%Y"))
  pas_unit<-(df_time_variable$time[2]-df_time_variable$time[1])

  # Dates correspondantes ---------------------------------------------------
 
  par(mfrow=c(3,4))
  liste_q_chosen<-rep(NA,37)
  liste_pu_chosen<-rep(NA,37)
  # Travail par temps pour les exts -----------------------------------------
  par(mfrow=c(1,1))

  # # Résultat avec arima -------------------------------------------
  # ######
  ## Juxtaposer les différentes séries.

  # Test avec un AR(1) en chaque temps --------------------------------------
  Obs_nouveau<-matrix(NA,nrow = nrow(obs),ncol = ncol(obs))
  Df_reg_ar<-matrix(NA,nrow=ncol(obs),ncol=P+1)
  Mat_ind_ext_p_valeur_Reiss<-list()
  Mat_ind_ext_p_valeur_Falk<-list()
  Mat_ind_ext_pval_Falk_brut<-list()
  Mat_Kendall<-list()
  Mat_borne<-rep(NA,ncol(obs))
  Mat_copula<-list()
  Mat_dep_measure<-list()
  Mat_copula_brut<-list()
  for(t in c(1:ncol(obs))){
    series<-obs[,t]
    modele<-arima(series,order=c(P,D,Q))
    Df_reg_ar[t,]<-unlist(modele$coef)
    residus<-modele$residuals
    Obs_nouveau[,t]<-residus
  
  }
  par(mfrow=c(1,2))
  for(ell in c(1)){
    series<-obs[,ell]
    ACF_alt<-acf(x=series,lag=400,plot=FALSE)
    plot(ACF_alt,main="",cex.lab=1.5,ylim=c(-1,1))
    PACF_alt<-pacf(x=series,lag=400,plot=FALSE)
    plot(PACF_alt,main="",cex.lab=1.5,ylim=c(-1,1))
  }
  par(mfrow=c(1,1))
  print("test")
  constante<-(1-Df_reg_ar[1,1])*Df_reg_ar[1,2]
  print(constante)
  print(c(obs[1,1],obs[2,1],Obs_nouveau[2,1]))
  print("fin test")
  
  Df_reg_ar<-as.data.frame(Df_reg_ar)
  Nom_cols_modele_lin<-c()
  for(j in c(1:P)){
    Nom_cols_modele_lin<-c(Nom_cols_modele_lin,paste0("coeff_ar",j))
  }
  colnames(Df_reg_ar)<-c(Nom_cols_modele_lin,"Intercept")
  nom_fichier<-paste0("differenciation_travail/",nom,"_Reg_AR",P,".csv")
  write.csv(file = nom_fichier,x=Df_reg_ar)
  
  graph_ar1<-melt(t(Df_reg_ar))
  colnames(graph_ar1)<-c("variable","temps","valeur")
  GG_graph_AR1<-ggplot(data=graph_ar1,aes(x=temps,y=valeur,col=variable))+
    facet_wrap(~variable,scales = "free_y")+
    geom_line()+
    geom_point()+
    ggtitle(paste0("Evolution du modèle AR(",P,") pour ",nom))+
    labs(col="Légende")+
    labs(caption=paste0("n=",nrow(Obs_nouveau)))
  # labs pour changer la légende. 
  print(GG_graph_AR1)
  
  # Analyse corrélation -----------------------------------------------------
  series_pos_AR<-c(t(Obs_nouveau))
  ACF_post_ar<-acf(series_pos_AR,lag=1000,plot=FALSE)
  PACF_post_ar<-pacf(series_pos_AR,lag=1000,plot=FALSE)
  
  series_complete<-as.vector(t(obs))
  titre<-paste0("Resultat obtenu avec (p=",P,", d=",D,", q=",Q,") en chaque temps")
  
  par(mfrow=c(1,2))
  # Analyse chaque temps ---------------------------------------------------
  vect_temps<-c(19)
  for(col in vect_temps){
    
    residus<-Obs_nouveau[,col]
    # travail avec ar en chaque temps -----------------------------------------
    ACF_alt<-acf(x=Obs_nouveau[,col],lag=1000,plot=FALSE)
    plot(ACF_alt,main="",cex.lab=1.5,
         ylim=c(-1,1))
    PACF_alt<-pacf(x=Obs_nouveau[,col],lag=1000,plot=FALSE)
    plot(PACF_alt,main="",cex.lab=1.5,
         ylim=c(-1,1))
    # Rajout caption pour détail nombre d'individus ds base. 
    #mtext(paste0("n=",nrow(Obs_nouveau)), side=1, line=3, adj = 1, cex=0.8)
    
    # mtext(paste0("Correlations for the residuals of ",nom_variable," at t=",col," with an AR(",P,")"), 
    #       line=-1,outer = TRUE,col = "red")
    # Test de Reiss-extRemes -----------------------------------------------------------
    # Test de Falk_POT --------------------------------------------------------
    
    Fct_test_ext_Reiss<-function(lag,series,C_th){
      fin<-length(series)-lag
      x_<-series[c(1:fin)]
      debut<-lag+1
      x_lag1<-series[c(debut:length(series))]
      #Reiss
      Test_lag1<-extRemes::taildep.test(x_,x_lag1,cthresh = C_th)
      return(Test_lag1)
    }
    # Mat_ind_ext_p_valeur_Reiss[[paste0("lag_",col)]]<-sapply(c(1:50),Fct_test_ext_Reiss,
    #                                                        series=residus,C_th=C_thchosen)
    # 
    # Fct_test_ext_Falk<-function(lag,series,C){
    #   fin<-length(series)-lag
    #   x_<-series[c(1:fin)]
    #   debut<-lag+1
    #   x_lag1<-series[c(debut:length(series))]
    #   #Falk
    #   Test_lag1<-POT::tailind.test(data = cbind(x_,x_lag1),
    #                                emp.trans =TRUE,c=C)$stats[,2]
    #   return(Test_lag1)
    # }
    # Fct_Kendall_coeff<-function(lag,series){
    #   fin<-length(series)-lag
    #   x_<-series[c(1:fin)]
    #   debut<-lag+1
    #   x_lag1<-series[c(debut:length(series))]
    #   # coefficient -------------------------------------------------------------
    #   Kendall<-cor.test(x_,x_lag1,method="kendall")
    #   return(Kendall)
    # }
    # Fct_analyse_copula<-function(lag,series,t){
    #   print(paste0("Copula analysis for lag=",lag))
    #   fin<-length(series)-lag
    #   x_<-series[c(1:fin)]
    #   debut<-lag+1
    #   x_lag1<-series[c(debut:length(series))]
    #   Unifs<-VineCopula::pobs(cbind(x_,x_lag1))
    #   Family_chosen<-VineCopula::BiCopSelect(Unifs[,1],Unifs[,2],indeptest = TRUE)
    #   return(Family_chosen)
    # }
    Fct_dep_measure<-function(lag,series,boot_bool,t){
      fin<-length(series)-lag
      x_<-series[c(1:fin)]
      debut<-lag+1
      x_lag1<-series[c(debut:length(series))]
      par(mfrow=c(1,2))
      resultat_chi<-POT::chimeas(data = cbind(x_,x_lag1),
                                 boot = boot_bool,
                                 ask = FALSE,which=1,cex.lab=1.5)
      abline(h=0,col="red")
      resultat_chi<-POT::chimeas(data = cbind(x_,x_lag1),
                                 boot = boot_bool,
                                 ask = FALSE,
                                 which=2,
                                 ylabs = rep(NA,2),
                                 cex.lab=1.5)
      abline(h=0,col="red")
      mtext(expression(bar(chi)),cex=1.5,las=2,
            outer=TRUE,line=-10.5)
      methode<-ifelse(test = boot_bool,yes = "boot", "delta")
      #title(paste0("Coefficients for a lag of ",lag, " for ", nom_variable," at (t=",t,") with the method ",methode),outer = TRUE,line = -2)
      par(mfrow=c(1,1))
      return(resultat_chi)
    }
    # 
    # Mat_ind_ext_p_valeur_Falk[[paste0("lag_",col)]]<-sapply(c(1:50),Fct_test_ext_Falk,
    #                                                       series=residus,C=C_Falk)
    # Mat_ind_ext_pval_Falk_brut[[paste0("lag_",col)]]<-sapply(c(1:50),Fct_test_ext_Falk,
    #                                                        series=series,C=C_Falk)
    # #Mat_copula[[paste0("lag_",col)]]<-sapply(c(1:10),Fct_analyse_copula,series=residus)
    # #Mat_copula_brut[[paste0("lag_",col)]]<-sapply(c(1:10),Fct_analyse_copula,series=series)
    # # Dependence measure ------------------------------------------------------
    Mat_dep_measure[[paste0("lag_",col)]]<-sapply(c(1:10),Fct_dep_measure,
                                                series=residus,boot_bool=FALSE,t=col)
    
  }
  par(mfrow=c(1,1))
  if(nom_variable=="S"){
    write.csv(x =Obs_nouveau,file = paste0("residus_clust/","Surcote","_residus.csv"))
  }
  else{
    write.csv(x =Obs_nouveau,file = paste0("residus_clust/",nom_variable,"_residus.csv"))
    
  }
  # Export du clustering ---------------------------------------------------
  #write.csv(x=obs_clust,file=paste0("residus_clust/clust_",nom_variable,".csv"))
  
  par(mfrow=c(1,2))
  
  # ### Normes
  NL2_deseason<-apply(Obs_nouveau,MARGIN = 1,FUN = calcul_norme_L2)
  
  plot(acf(NL2_deseason,lag=1000,plot=FALSE),main="ACF",ylim=c(-1,1))
  plot(pacf(NL2_deseason,lag=1000,plot=FALSE),main="PACF",ylim=c(-1,1))
  mtext(paste0("n=",nrow(Obs_nouveau)), side=1, line=3, adj = 1, cex=0.8)
  mtext(paste0("Corrélations pour la norme L2 de ",nom_variable," avec un AR(",P,")"), 
        line=-1,outer = TRUE,col = "red")

  # dépendance asymptotique -------------------------------------------------
  for( lag_norme in c(1:10)){
    par(mfrow=c(1,2))
    fin<-length(NL2_deseason)-lag_norme
    serie_1<-NL2_deseason[1:fin]
    series_lag<-NL2_deseason[lag_norme+1:length(NL2_deseason)]
    resultat_chi<-POT::chimeas(data = cbind(serie_1,series_lag),
                               boot = FALSE,ask=FALSE,which=1)
    abline(h=0,col="red")
    
    resultat_chi<-POT::chimeas(data = cbind(serie_1,series_lag),
                               boot = FALSE,ask=FALSE,which=2)
    abline(h=0,col="red")
    methode<-ifelse(test = FALSE,yes = "boot", "delta")
    title(paste0("Coefficients pour un décalage de ",lag_norme, " pour la norme de ", nom_variable," via ",methode),
          outer = TRUE,line = -2)
    par(mfrow=c(1,1))
  }
  
  plot(acf(VECTEUR_normeL2,lag=1000,plot=FALSE),main=paste0("ACF NL2 orig (",nom_variable,")"),ylim=c(-1,1))
  plot(acf(NL2_deseason,lag=1000,plot=FALSE),main=paste0("ACF NL2 résidus de (",nom_variable,")"),ylim=c(-1,1))
  mtext(paste0("n=",nrow(Obs_nouveau)), side=1, line=3, adj = 1, cex=0.8)
  mtext(titre, line=-1,outer = TRUE,col = "red")
  
  plot(pacf(VECTEUR_normeL2,lag=1000,plot=FALSE),main=paste0("PACF NL2 orig (",nom_variable,")"),ylim=c(-1,1))
  plot(pacf(NL2_deseason,lag=1000,plot=FALSE),main=paste0("PACF NL2 résidus de (",nom_variable,")"),ylim=c(-1,1))
  mtext(paste0("n=",nrow(Obs_nouveau)), side=1, line=3, adj = 1, cex=0.8)
  mtext(titre, line=-1,outer = TRUE,col = "red")
  
  plot(acf(NL2_deseason,lag=1000,plot=FALSE),main=paste0("ACF NL2 résidus de (",nom_variable,")"))
  plot(pacf(NL2_deseason,lag=1000,plot=FALSE),main=paste0("PACF NL2 résidus de (",nom_variable,")"))
  mtext(paste0("n=",nrow(Obs_nouveau)), side=1, line=3, adj = 1, cex=0.8)
  mtext(titre, line=-1,outer = TRUE,col = "red")
  
  return(list("df"=df_time_variable,"df_clust"=NULL,
              "dates_ext"=NULL,"clusters_id"=NULL,
              "fonction_seuil"=liste_q_chosen,"df_deason"=Obs_nouveau,
              "Matrice_test_Reiss"=Mat_ind_ext_p_valeur_Reiss,
              "Matrice_test_Falk"=Mat_ind_ext_p_valeur_Falk,
              "Matrice_test_Falk_brut"=Mat_ind_ext_pval_Falk_brut,
              "Matrice_chi"=Mat_dep_measure,
              "Matrice_Copula"=Mat_copula, 
              "Matrice_Copula_brut"=Mat_copula_brut))
}

nom<-"Surcote"
vecteur<-list()
S<-37
C_TH<--0.5
CF<--0.1
for(p_sar in c(1)){
  R<-fonction_arima_variable(coeurs = coeurs,nom_variable = nom,opt_tend = TRUE,type_donnees = "HIVER",P=p_sar,
                             Q=0,D=0,p=0,q=0,S=S, pas_clust=1,C_thchosen=C_TH,C_Falk=CF)
}

# analyse des dépendances ou non presentes --------------------------------

# Analyse des dépendances -------------------------------------------------
Chi<-R$Matrice_chi
Chi_temps_1<-Chi$lag_1
elt_lag11<-Chi_temps_1[,1]$chi
plot(elt_lag11[1,],elt_lag11[2,],type="l")
lines(elt_lag11[1,],elt_lag11[3,],type="l",col="red")
lines(elt_lag11[1,],elt_lag11[4,],type="l",col="red")

# Données brutes ----------------------------------------------------------
lag1_brut<-R$Matrice_test_Falk_brut$lag_1
lag1_brut
# Données résidus ---------------------------------------------------------

test_Reiss<-R$Matrice_test_Reiss
test_Falk<-R$Matrice_test_Falk
test_Falk$lag_1
names(test_Reiss)
vect_repere<-c("1","19","37")

pval<-matrix(NA,nrow=length(vect_repere),ncol=ncol(test_Reiss$lag_1))
pval_Falk<-matrix(NA,nrow=length(vect_repere)*4,ncol=ncol(test_Reiss$lag_1))
Taux<-matrix(NA,nrow=length(vect_repere),ncol=ncol(test_Reiss$lag_1))

for(j in c(1:length(vect_repere))){
  elt<-vect_repere[j]
  nom_col<-paste0("lag_",elt)
  mat_lag<-test_Reiss[[nom_col]]
  pval[j,]<-as.numeric(unlist(mat_lag["p.value",]))
  falk_j<-test_Falk[[nom_col]]
  debut<-(j-1)*nrow(falk_j)+1
  fin<-j*nrow(falk_j)
  pval_Falk[c(debut:fin),]<-falk_j
  Taux[j,]<-matrix(unlist(mat_lag["parameter",]),ncol=4,byrow = TRUE)[,4]
}

# Travail sur le test de Falk ---------------------------------------------
###############

resultat_melteur<-melt(t(pval_Falk))
resultat_melteur$value
L<-length(rownames(falk_j))
temps<-sort(rep(vect_repere,L))
type_test<-rep(sapply(X=rownames(falk_j),FUN = function(x){return(rep(x,50))}),3)

colnames(resultat_melteur)<-c("lag","temps","valeur")
resultat_melteur$temps<-as.character(temps[resultat_melteur$temps])
resultat_melteur$test<-type_test
resultat_melteur$hyp_rejetee<-as.character(ifelse(resultat_melteur$valeur<0.05,0,1))

resultat1<-resultat_melteur[which(resultat_melteur$test=="KS"),]
ggplot(data=resultat1,aes(x=lag,y=valeur,col=test))+
  facet_wrap(~temps)+
  geom_point(aes(shape=hyp_rejetee))+
  geom_hline(aes(yintercept = 0.05,col="seuil_confiance"))+
  ggtitle(paste0("p valeur de l'hypothèse d'indépendance extrémale (Falk) pour ",nom))+
  labs(col="Légende",caption=paste0("valeur de c=",CF))

resultat1<-resultat_melteur[which(resultat_melteur$test=="ChiSq"),]
ggplot(data=resultat1,aes(x=lag,y=valeur,col=test,pch=hyp_rejetee))+
  facet_wrap(~temps)+
  geom_point()+
  geom_hline(aes(yintercept = 0.05,col="seuil_confiance"))+
  ggtitle(paste0("p valeur de l'hypothèse d'indépendance extrémale (Falk) pour ",nom))+
  labs(col="Légende",caption=paste0("valeur de c=",CF))

resultat1<-resultat_melteur[which(resultat_melteur$test=="Fish"),]
ggplot(data=resultat1,aes(x=lag,y=valeur,col=test,pch=hyp_rejetee))+
  facet_wrap(~temps)+
  geom_point()+
  geom_hline(aes(yintercept = 0.05,col="seuil_confiance"))+
  ggtitle(paste0("p valeur de l'hypothèse d'indépendance extrémale (Falk) pour ",nom))+
  labs(col="Légende",caption=paste0("valeur de c=",CF))

resultat1<-resultat_melteur[which(resultat_melteur$test=="NP"),]
ggplot(data=resultat1,aes(x=lag,y=valeur,col=test,pch=hyp_rejetee))+
  facet_wrap(~temps)+
  geom_point()+
  geom_hline(aes(yintercept = 0.05,col="seuil_confiance"))+
  ggtitle(paste0("p valeur de l'hypothèse d'indépendance extrémale (Falk) pour ",nom))+
  labs(col="Légende",caption=paste0("valeur de c=",CF))

resultat1<-resultat_melteur[which(resultat_melteur$test=="KS"),]
ggplot(data=resultat1,aes(x=lag,y=valeur,col=test,pch=hyp_rejetee))+
  facet_wrap(~temps)+
  geom_point()+
  geom_hline(aes(yintercept = 0.05,col="seuil_confiance"))+
  ggtitle(paste0("p valeur de l'hypothèse d'indépendance extrémale (Falk) pour ",nom))+
  labs(col="Légende",caption=paste0("valeur de c=",CF))


# Retour sur le test de Reiss --------------------------------------------
#####################
df_proportion<-melt(t(Taux))
colnames(df_proportion)<-c("lag","temps","proportion")
df_proportion$temps<-vect_repere[df_proportion$temps]

df<-cbind.data.frame(melt(t(pval)))
colnames(df)<-c("lag","temps","p.valeur")
df$temps<-vect_repere[df$temps]

ggplot(data=df,aes(x=lag,y=p.valeur,col=temps))+
  facet_wrap(~temps)+
  geom_point()+
  geom_hline(aes(yintercept = 0.01,col="seuil_confiance"))+
  ggtitle(paste0("p valeurs du test de Reiss pour ",nom))+
  labs(caption=paste0("cth=",C_TH))
  

ggplot(data=df_proportion,aes(x=lag,y=proportion,col=temps))+
  facet_wrap(~temps)+
  geom_point()+
  geom_hline(aes(yintercept = 10,col="borne_conseille"))+
  geom_hline(aes(yintercept = 15,col="borne_conseille"))+
  ggtitle(paste0("Proportion d'excédents du test de Reiss pour ",nom))+
  labs(caption=paste0("cth=",C_TH))
  

par(mfrow=c(1,1))
resultat<-c(t(R$df_deason))
which(R$dates_ext=="17/12/1989")
resultat_nb_exts<-table(substr(R$dates_ext,7,10))
plot(resultat_nb_exts)

time<-1
###### Copula analysis of dependence.

pval_raw<-c()
emptau_raw<-c()
Base_travail<-R$Matrice_Copula_brut[[paste0("lag_",time)]]
for(j in c(1:ncol(Base_travail))){
  emptau_raw<-c(emptau_raw,Base_travail[,j]$emptau)
  pval_raw<-c(pval_raw,Base_travail[,j]$p.value.indeptest)
}
df_pval_raw<-cbind.data.frame(c(1:length(pval)), 
                          pval_raw, emptau_raw)
colnames(df_pval_raw)<-c("h", "pval","kendall")
threshold_conf<-0.05
ggplot(data=df_pval_raw,aes(x=h,y=pval))+
  geom_point()+
  geom_hline(aes(yintercept =threshold_conf,col="threshold"))+
  annotate("text", x = 0, y = 2*threshold_conf, 
           label=threshold_conf)+
  xlab("Lag h")+ylab("p value")+
  labs(color="Legend")+
  ggtitle(paste0("p value of the independence test at t=",time, " for ",nom))

pval<-c()
emptau<-c()
###### Copula analysis of dependence_residuals
Base_travail<-R$Matrice_Copula[[paste0("lag_",time)]]
for(j in c(1:ncol(Base_travail))){
  emptau<-c(emptau,Base_travail[,j]$emptau)
  pval<-c(pval,Base_travail[,j]$p.value.indeptest)
}
df_pval<-cbind.data.frame(c(1:length(pval)), 
                          pval,emptau)
colnames(df_pval)<-c("h", "pval","kendall")
ggplot(data=df_pval,aes(x=h,y=pval))+
  geom_point()+
  geom_hline(aes(yintercept =threshold_conf,col="threshold"))+
  annotate("text", x = 0, y = 2*threshold_conf, 
           label=threshold_conf)+
  xlab("Lag h")+ylab("p value")+
  labs(color="Legend")+
  ggtitle(paste0("p value of the independence test at t=",time, " for ",nom, " (residuals)"))

ggplot(data=df_pval,aes(x=h,y=kendall,col="residuals"))+
  geom_point(pch=2)+
  geom_point(data = df_pval_raw,aes(x=h,y=kendall,col="data"))+
  xlab("Lag h")+ylab("Kendall coefficient")+
  labs(color="Legend")+
  geom_hline(yintercept = 0)+
  ggtitle(paste0("Kendall coefficient at t=",time," for ", nom))
df_pval_raw 
# Residus modèle ----------------------------------------------------------
Series_R<-as.vector(t(R$df_deason))
plot(density(Series_R))
Res_pos<-subset(Series_R,Series_R>3)
plot(density(Res_pos))

# t<-1
# par(mfrow=c(2,2))
# analyse_post_choix<-function(data,t,f_seuil){
#   obs_temps<-as.data.frame(list("obs"=data[,t],"time"=1:nrow(data)))
#   Q1<-quantile(data[,t],0.1)
#   Q2<-quantile(data[,t],0.95)
#   title(paste0("Graphiques pour t=",t),outer=TRUE,line=-3)
#   
#   valeurs<-POT::tcplot(data[,t],ask = FALSE,u.range = c(Q1,Q2),which=1)
#   abline(v=f_seuil[t],col="green")
#   valeurs<-POT::tcplot(data[,t],ask = FALSE,u.range = c(Q1,Q2),which=2)
#   abline(v=f_seuil[t],col="green")
#   POT::mrlplot(data = data[,t],u.range = c(Q1,Q2))
#   abline(v=f_seuil[t],col="green")
#   #POT::diplot(data = obs_temps,nt=100)
#   #abline(v=f_seuil[t],col="green")
#   
#   ### Analyse gpd #####
#   observations_GPD_question<-data[which(data[,t]>f_seuil[t]),t]
#   
# }
# resultat_temps_choisi<-analyse_post_choix(data=R$df_clust,t=t,f_seuil=R$fonction_seuil)

# Fin code ----------------------------------------------------------------
stopCluster(coeurs)

