rm(list=ls())
set.seed(133)
par(mfrow=c(1,1))

#### Import des packages. #####
require(dplyr)
require(ggplot2)
require(reshape2)
require(FactoMineR)
require(fda)
require(parallel)
require(VineCopula)
require(mvPot)
require(SpatialExtremes)
### Source des fonctions. ####
source("fonctions_Pareto.R")
source("Modele_ss_tend/Travail_angles.R")
source("fonctions.r")
source("fonctions_r_pareto_EVOLTAU.r")
F_names<-ls()
nb_coeurs<-detectCores()-4


# Ouverture des clusters --------------------------------------------------
coeurs<-makeCluster(nb_coeurs)
clusterExport(coeurs,varlist=F_names)
# (I) Resultat des simulations ----------------------------------------------------------------
# #### --------------------------------------------------------------------

nom_variable<-"Hs"
TS_GPD<-"POT"
typeD<-"HIVER"
BY<-c(0,13)

fonction_analyse_gamma<-function(coeurs,p_norm,nom,n.dens,opt_tend,Tau,
                                 BORNES_Y=c(0,30),type_donnees){

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
  plot(colMeans(obs),main="Fonction moyenne")
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
  Tendance_ounon<-ifelse(opt_tend,"On enlève la tendance","On conserve la tendance")
  Vecteur_DIM<-1:ncol(obs)

  # loi suivie par la norme L2 ----------------------------------------------
  L2_norme<-apply(X = obs,MARGIN = 1,FUN = calcul_norme_L2)
  seuil_L2_orig<-quantile(L2_norme,1-p_norm)
  NPY<-nrow(obs)/37
  modele_ev<-extRemes::fevd(x = L2_norme,
                            threshold = seuil_L2_orig,type="GP",
                            time.units = paste0(NPY,"/year"))
  print(modele_ev$results$par)
  L2_exts<-subset(L2_norme,L2_norme>seuil_L2_orig)
  print(NPY)
  plot(modele_ev)
  periods_years<-c(2,5,10,20,50,80,
                   100,120,200,
                   250,300,500,
                   800)
  RL_plot<-plot(modele_ev,type = "rl",
       main=paste0("Niveau de retour par EVA pour la norme de ",nom),
       rperiods=c(periods_years))
  empirique<-RL_plot$empirical
  Modele_vs_emp<-RL_plot$model
  Base<-cbind.data.frame(Modele_vs_emp[,1],Modele_vs_emp[,2],Modele_vs_emp[,3])
  colnames(Base)<-c("borne_inf","estimateur","borne_sup")
  print(nom)
  last_letter<-substr(nom,length(nom),length(nom))
  print(length(nom))
  last_letter_stuff<-ifelse(last_letter=="s","'","'s")
  Cols_<-c("data"="blue","GPD"="green","confidence_bands"="darkblue")
  GG_orig<-ggplot(data = Base,aes(x=periods_years,y=estimateur,col="GPD"))+
    geom_line()+
    geom_ribbon(data =Base,mapping = aes(ymin=borne_inf,ymax=borne_sup,col="confidence_bands"),alpha=0.15,fill="grey",
                linetype = "dashed")+
    geom_point(data=empirique,aes(x=transformed.period,
                                  y=sorted.level,col="data"))+
    ylab("L2 norm")+
    xlab("Period (years)")+
    scale_x_continuous(transform = "log")+
    scale_color_manual(values=Cols_)+
    labs(colour="Legende",caption=paste0("Number of extremes=",length(L2_exts),", npy=",round(NPY),
                                         ", threshold=",round(seuil_L2_orig,2),
         ", total number of observations=",length(L2_norme)))+
    ggtitle(paste0("Return level plot for ",nom,last_letter_stuff, " norm"))
  print(GG_orig)
    
  # Comparaison marginales Pareto, modele de melange ------------------------
  
  Resultat_P<-as.data.frame(t(sapply(Vecteur_DIM,f_marginales_all_Pareto,donnees=obs,p_u=p_u,n.dens=n.dens,
                                        donnees_a_transformer=obs)))
  K<-Resultat_P$kernel_dens
  PARETO<-Resultat_P$obs
  vecteur_pareto<-do.call(cbind.data.frame,PARETO)
  NL2<-apply(X = vecteur_pareto,MARGIN = 1,FUN = calcul_norme_L2)
  ul<-quantile(NL2,probs=1-p_norm)
  vecteur_pareto_conserves<-subset(vecteur_pareto,NL2>ul)
  Pareto_NL2<-subset(NL2,NL2>ul)
  print(length(NL2))
  scale<-ul
  shape<-1
  Gamma_ML<-shape
  seuil<-ul
  p_normale<-sort(seq.int(from = 0.01,to = 0.50,length.out = 30)^(-1))
  fonction_GPD_exp_quantile<-function(scale,shape,seuil,p){
    part1<-p^(-shape)-1
    return((part1/shape)*scale+seuil)
  }
  Ninter<-1000
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
  simulations_for_int<-QGPDpareto(ppoints(n =Ninter))
  valeurs_associees_theorique<-sapply(p_normale^(-1),FUN = fonction_GPD_exp_quantile,scale=scale,shape=Gamma_ML,seuil=seuil)
  valeurs_associees_empirique<-as.numeric(quantile(Pareto_NL2,1-p_normale^(-1)))
  Quantiles_full_theor<-extRemes::qqplot(simulations_for_int,simulations_for_int,make.plot=FALSE)
  upper<-Quantiles_full_theor$qdata$upper
  lower<-Quantiles_full_theor$qdata$lower
  
  indices<-which((is.na(lower)==FALSE)&(is.na(upper)==FALSE))
  B_supKs<-as.numeric(quantile(upper[indices],1-p_normale^(-1)))
  B_infKs<-as.numeric(quantile(lower[indices],1-p_normale^(-1)))
  aspect_extreme<-cbind.data.frame(valeurs_associees_empirique,valeurs_associees_theorique)
  colnames(aspect_extreme)<-c("empirique","GPD")
  Melteur<-melt(t(aspect_extreme))
  colnames(Melteur)<-c("source","probabilite_depassement","valeur")
  Melteur$sup_KS<-B_supKs[Melteur$probabilite_depassement]
  Melteur$inf_KS<-B_infKs[Melteur$probabilite_depassement]
  Melteur$probabilite_depassement<-p_normale[Melteur$probabilite_depassement]
  Cols_<-c("empirique"="blue","GPD"="green","Percentile_5_95_KS"="darkblue")
  GGreturn_log<-ggplot(data=Melteur,aes(x =probabilite_depassement,y=valeur,color=source,group=interaction(source)))+
    geom_point()+
    geom_line()+
    scale_x_continuous(trans = "log10")+
    labs(caption=paste0("simulations_GPD_KS=",Ninter,", gamma=", Gamma_ML,", seuil=", round(ul,2),", scale=",round(ul,2)))+
    scale_color_manual(values=Cols_)+
    xlab("1/p")+
    geom_ribbon(data =Melteur,mapping = aes(ymin=inf_KS,ymax=sup_KS,col="Percentile_5_95_KS"),alpha=0.15,fill="grey", linetype = "dashed")+
    ggtitle(paste0("Comportement de la norme L2 de T(",nom_variable,") (repère log(x)) pour les excédents"))
  print(GGreturn_log)
  
  Matrice_couples<-list()
  L<-ncol(vecteur_pareto)
  z<-1
  vecteur_distances<-c()
  for(j in 1:L){
    for (i in 1:L){
      Matrice_couples[[z]]<-c(i,j)
      vecteur_distances<-c(vecteur_distances,abs(i-j))
      z<-z+1
    }
  }
  Inv_Extremo_emp<-inv_extremogram_empirique(Matrice_couples=Matrice_couples,individus_select = vecteur_pareto_conserves,Tau = 1/(p_u))
  return(list("extremo"=Inv_Extremo_emp,"distance"=vecteur_distances,"individus"=vecteur_pareto_conserves))
}

lancement_code_Gamma_analyse<-function(coeurs,p_norm,nom,n.dens,opt_tend,type_seuils,BORNES_Y=c(0,30),type_donnees){
  e_<-read.table(paste0("SEUILS_THRESHR_",nom_variable,"_",type_donnees,"_",type_seuils,".txt"))
  Seuils_<-e_$x
  plot(c(1:37),Seuils_,main=paste0("Seuils utilisés pour ",nom_variable," en ",tolower(type_donnees)),xlab="Temps",ylab="Seuil",type="o")
  plot(c(1:37),cummean(Seuils_),main=paste0("Seuils utilisés pour ",nom_variable," en ",tolower(type_donnees)," (moyenne cumulée)"),xlab="Temps",ylab="Seuil",type="o")
  resultat<-fonction_analyse_gamma(coeurs=coeurs,p_norm=p_norm,nom=nom,n.dens=n.dens,opt_tend=opt_tend,Tau=cummean(Seuils_),BORNES_Y=BORNES_Y,type_donnees=type_donnees)
  stopCluster(coeurs)
  return(resultat)
}

extremo<-lancement_code_Gamma_analyse(coeurs = coeurs,p_norm = 0.05,nom = nom_variable,n.dens = 2^(13),opt_tend = TRUE,type_seuils=TS_GPD,BORNES_Y = BY,type_donnees = typeD)
Inds_extremes<-extremo$individus

NL2_extremes<-apply(X = Inds_extremes,MARGIN = 1,FUN = calcul_norme_L2)
Angle_extremes<-t(t(Inds_extremes)%*%diag(NL2_extremes^(-1)))
extremo_sous<-extremo$extremo
df<-cbind.data.frame(extremo$extremo,extremo$distance)

colnames(df)<-c("valeur","delta")
require(dplyr)

resultat_delta<-df %>% group_by(delta) %>% summarise(avg=mean(valeur))
Avg<-resultat_delta$avg

# Brown_Resnick -----------------------------------------------------------
#####################

# semi_vario<-2*(qnorm(1-(Avg/2)))^(2)
# plot(resultat_delta$delta/37,semi_vario,type="o",xlab="Ecart",ylab="Valeur",main=paste0("Approximation du semi-variogramme pour ",nom_variable))
# summary(resultat_delta)
# Modèle semi-variogramme à essayer.  ----------------------------------------------

# Géométrique -------------------------------------------------------------
####
sigma_moins_cov<-(qnorm(Avg/2))^(2)*2
print(length(sigma_moins_cov))
# En 0, vient de l'égalité entre c(h) et sigma carré.
NC<-"whitmat"
fonction_cout<-function(parametre,cible,vecteur_delta){
  if(parametre[1]<0|parametre[2]<0|parametre[3]<0|parametre[3]>2){return(1e50)}
  sigma<-parametre[1]
  lambda<-parametre[2]
  kappa<-parametre[3]
  f_cov<-function(parametre,nom_cov,h){
    f_ch<-covariance(nugget = 0, sill = sigma, range =lambda, smooth = kappa, cov.mod =nom_cov,plot = FALSE)
    return(f_ch(h))
  }
  ch<-sapply(X=vecteur_delta,FUN=f_cov,parametre=parametre,nom_cov=NC)
  return(mean((sigma-ch-cible)^2))
}
Params_nuls<-c(1,1,0.5)
D<-resultat_delta$delta/nrow(resultat_delta)
Optimisation_Geometric<-optim(par = Params_nuls,fn = fonction_cout,cible=sigma_moins_cov,vecteur_delta=D,method="Nelder-Mead",control=list(maxit = 1000))
Parametres<-Optimisation_Geometric$par
f_optim<-SpatialExtremes::covariance(nugget = 0, sill = Parametres[1], range =Parametres[2], smooth = Parametres[3], cov.mod =NC,plot = FALSE)
retour<-sapply(X = D,FUN = f_optim)

plot(D,sigma_moins_cov,main="Résultat modélisation Geometric (théorique en rouge)")
lines(D,Parametres[1]-retour,col="red")
print(Parametres)
# Simulation_GP -----------------------------------------------------------
Realisations<-SpatialExtremes::rgp(n = 2500,coord = c(1:37)/37,mean = 0,sill = Parametres[1],cov.mod = "whitmat",nugget = 0,range = Parametres[2],smooth = Parametres[3],grid = FALSE)
Geometric<-exp(Realisations-Parametres[1]/2)
NL2<-apply(X = Geometric,MARGIN = 1,FUN = calcul_norme_L2)
realisations<-subset(Geometric,NL2>1)
NL2_ext<-subset(NL2,NL2>1)
Angle<-t(t(realisations)%*%diag(NL2_ext^(-1)))
print(nrow(Angle))
par(mfrow=c(1,2))
matplot(t(Angle),type="l",ylim=c(0,2.55))
matplot(t(Angle_extremes),type="l",ylim=c(0,2.55))

# Simulation d'une Pareto -------------------------------------------------


