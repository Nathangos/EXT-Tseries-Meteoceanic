rm(list=ls())
set.seed(133)

### Source des fonctions. ####
source("fonctions/fonctions_r_pareto_EVOLTAU.R")

#### Packages. #####

require(patchwork)
require(corrplot)
require(ggplot2)
require(dplyr)
require(reshape2)
require(FactoMineR)
require(fda)
require(extRemes)
require(parallel)
require(POT)

source("fonctions/fonctions.R")
F_names<-ls()
nb_coeurs<-detectCores()-4
coeurs<-makeCluster(nb_coeurs)
clusterExport(coeurs,varlist=F_names)
periode<-c(1979,2016)
annees<-seq.int(periode[1],periode[2],by=1)
dataFiles<-c()

for(elt in annees){
  commentaire<-paste0("donnees_transformees/evt_",elt,"*.csv")
  file_annee <- parLapply(cl=coeurs,Sys.glob(commentaire), read.csv)
  dataFiles<-c(dataFiles,file_annee)
}


# Extraction Johanna ------------------------------------------------------
dates_J<-c("08/03/2008","09/03/2008","10/03/2008")
Analyse_univs<-function(df,nom_variable,coeurs,Nb_sampling,period_analysis){

  # Extract time ------------------------------------------------------------
  TEMPS<-parLapply(cl=coeurs,df,f_nom_variable,"Date")
  intermed<-t(do.call(cbind.data.frame,TEMPS))
  rownames(intermed)<-1:nrow(intermed)
  
  # Date pic de la marée.  --------------------------------------------------
  date_compar<-intermed[,19]
  colonne<-as.numeric(substr(date_compar,4,5))
  
  obs<-parLapply(cl=coeurs,df,f_nom_variable,nom_variable)
  matrice_<-matrix(NA,nrow=length(obs),ncol=37)
  for (col in c(1:37)){
    variable<-sapply(obs,FUN = function(t,serie){return(serie[t])},t=col)
    matrice_[,col]<-variable
  }
  Fonction_Btrap_entree_fonct<-function(nb_echs,size_ech,fonction_to_apply,data_to_analyse){
    fonction_by_ech<-function(data_to_analyse,size_ech,fonction_to_apply){
      reech<-sample(c(1:nrow(data_to_analyse)),size = size_ech,replace = TRUE)
      sub_ech<-data_to_analyse[reech,]
      return(apply(X = sub_ech,FUN = fonction_to_apply,MARGIN=2))
    }
    result<-replicate(nb_echs,expr =fonction_by_ech(data_to_analyse,size_ech,fonction_to_apply))
    return(result)
  }

  # Month_evol_variable_peak ------------------------------------------------------
  Nom_graph<-ifelse(test = nom_variable=="Surcote",yes = "S",
                    no=nom_variable)
  df_month_peak<-cbind.data.frame(matrice_[,19],colonne)
  colnames(df_month_peak)<-c("variable","month_o")
  month_winter<-c(1:3,9:12)
  df_month_peak$month_winter<-ifelse(df_month_peak$month_o%in%month_winter,
                                     yes ="Winter" ,no="Other")
  
  Nom_graph<-ifelse(nom_variable=="Surcote",yes = "S",
                    no=nom_variable)
  GG_violin<-ggplot(data=df_month_peak,aes(x=factor(month_o),y=variable))+geom_violin()+
    geom_boxplot(alpha=0.15,col="black")+
    ggtitle(paste0("Evolution of the distribution of ",Nom_graph))
  
  FQ25<-function(x){return(quantile(x,0.25))}
  FQ75<-function(x){return(quantile(x,0.75))}

  # M out of n bootstrap ----------------------------------------------------
  Mean_t<-Fonction_Btrap_entree_fonct(nb_echs=1000,size_ech=1000,
                                          fonction_to_apply = mean,data_to_analyse = matrice_)
  borne_sup_mean<-apply(X = Mean_t,MARGIN = 1,FUN = function(x){return(quantile(x,0.975))})
  borne_inf_mean<-apply(X = Mean_t,MARGIN = 1,FUN = function(x){return(quantile(x,0.025))})
  Mean_whole<-cbind.data.frame(borne_sup_mean,borne_inf_mean)
  colnames(Mean_whole)<-c("ymax","ymin")
  Mean_whole$variable<-rep("mean",nrow(Mean_whole))
  print(nrow(matrice_))
  # Quantile ----------------------------------------------------------------
  Quantile_t_25<-Fonction_Btrap_entree_fonct(nb_echs=1000,size_ech=1000,
                                          fonction_to_apply = FQ25,data_to_analyse = matrice_)
  borne_sup_Q25<-apply(X = Quantile_t_25,MARGIN = 1,FUN = function(x){return(quantile(x,0.975))})
  borne_inf_Q25<-apply(X = Quantile_t_25,MARGIN = 1,FUN = function(x){return(quantile(x,0.025))})
  Q25_whole<-cbind.data.frame(borne_sup_Q25,
                              borne_inf_Q25)
  colnames(Q25_whole)<-c("ymax","ymin")
  Q25_whole$variable<-rep("Q_25",nrow(Mean_whole))
  
  Quantile_t_75<-Fonction_Btrap_entree_fonct(nb_echs=1000,size_ech=1000,
                                             fonction_to_apply = FQ75,data_to_analyse = matrice_)
  borne_inf_Q75<-apply(X = Quantile_t_75,MARGIN = 1,FUN = function(x){return(quantile(x,0.025))})
  borne_sup_Q75<-apply(X = Quantile_t_75,MARGIN = 1,FUN = function(x){return(quantile(x,0.975))})
  Q75_whole<-cbind.data.frame(borne_sup_Q75,
                              borne_inf_Q75)
  colnames(Q75_whole)<-c("ymax","ymin")
  Q75_whole$variable<-rep("Q_75",nrow(Mean_whole))
  
  # standard deviation ------------------------------------------------------
  Sigma_t<-Fonction_Btrap_entree_fonct(nb_echs=1000,size_ech=1000,
                                       fonction_to_apply = sd,data_to_analyse = matrice_)
  sigma_sup<-apply(X = Sigma_t,MARGIN = 1,FUN = function(x){return(quantile(x,0.975))})
  sigma_inf<-apply(X = Sigma_t,MARGIN = 1,FUN = function(x){return(quantile(x,0.025))})
  Sigma_whole<-cbind.data.frame( sigma_sup,
                                 sigma_inf)
  colnames(Sigma_whole)<-c("ymax","ymin")
  Sigma_whole$variable<-rep("std",nrow(Mean_whole))
  
  # Stats_univariees -----------------------------------------------------------------
  Mean_whole$value<-apply(X = matrice_,MARGIN = 2,FUN = mean)
  Sigma_whole$value<-apply(X=matrice_,MARGIN = 2,FUN=sd)
  Q25_whole$value<-apply(X=matrice_,MARGIN = 2,FUN=function(x){return(quantile(x,0.25))})
  Q75_whole$value<-apply(X=matrice_,MARGIN = 2,FUN=function(x){return(quantile(x,0.75))})
  All_stats<-rbind.data.frame(Mean_whole,Q25_whole,Q75_whole,
                              Sigma_whole)
  All_stats$Time<-rep(c(1:37),4)
  All_stats$Time<-(All_stats$Time-19)*(1/6)
  GG_univ<-ggplot(All_stats,aes(x=Time,y=value,interaction=variable,col=variable))+
    geom_line()+
    facet_wrap(~variable,scales = "free")+
    #ggtitle(paste0("Univariate analysis per time for ",ifelse(nom_variable=="Surcote","S",
     #                                                         nom_variable)," (", period_analysis[1],",",period_analysis[2],")"))+
    xlab("Time (hour) with respect to tidal peak")+
    ylab("value (m)")+
    geom_ribbon(aes(ymin=ymin,ymax=ymax, 
                    col="confidence_band"),alpha=0.15,
                fill="grey", linetype = "dashed")+
    scale_color_manual(values=c("confidence_band"="darkblue",
                                "mean"="orange",
                                "std"="red",
                                "Q_25"="green",
                                "Q_75"="blue"))+
    labs(col="Legend")+
    theme(axis.title=element_text(size=15))
    
  print(GG_univ)

  # Extreme value -----------------------------------------------------------
  #TITRE<-paste0("Shape parameter of the L2 norm of ",
  #       Nom_graph)
  vect_k<-c(20:600)
  L2_raw<-apply(X = matrice_,MARGIN = 1,FUN = calcul_norme_L2)
  Graphics_estimators_gamma(series =L2_raw,Title_graphic =NULL ,
                               vect_k = vect_k,
                               NB_years = 37)

  # Tirage aléatoire --------------------------------------------------------
  tirage_aleatoire<-sample(c(1:nrow(matrice_)),size = Nb_sampling,replace = FALSE)
  SS_ech<-matrice_[tirage_aleatoire,]
  melteur_sample<-melt(t(SS_ech))
  colnames(melteur_sample)<-c("time","individual","value")
  melteur_sample$individual<-as.character(melteur_sample$individual)
  GG_sample<-ggplot(data = melteur_sample,aes(x=time,y=value,interaction=individual,col=individual))+
    geom_line()+
    guides(col="none")+
    xlab("Time")+
    ylab("Surge (m)")+
    theme(axis.title=element_text(size=15),
          legend.text=element_text(size=10))
  print(GG_sample)
  return(list("month"=df_month_peak,"MAT_for_deriv"=matrice_))
}
Nom_v<-"Surcote"
df<-Analyse_univs(df=dataFiles,coeurs=coeurs,nom_variable =Nom_v , 
                  Nb_sampling=20,period_analysis=periode)
# apply(df,MARGIN = 2,FUN = function(x){return(quantile(x,0.05))})
# # df
# # df_pred<-matrix(NA,nrow = nrow(df),
# #                 ncol=ncol(df))
# # LM_vect<-function(vect_y){
# #   mean_y<-mean(vect_y)
# #   time<-c(1:37)/36
# #   alpha<-mean((vect_y-mean_y)*(time-mean(time)))/var(time)
# #   beta<-mean_y-alpha*mean(time)
# #   L2<-calcul_norme_L2(abs(alpha*time+beta))
# #   return(L2)
# # }
# # Predictions_lin<-apply(df,MARGIN = 1,FUN = LM_vect)
# # Value<-apply(X = abs(df),MARGIN = 1,FUN = calcul_norme_L2)
# # summary(Value)
# # summary(Predictions_lin)
# # Inds<-which(Value!=0)
# # R_sq<-Predictions_lin[Inds]/Value[Inds]
# # 
# # 
# # # Coeff_<-(df[,37]-df[,1])/(36/37)
# # # Pred<-sapply(X = c(0:36),FUN = function(x){return(Coeff_*x)})/37+df[,1]
# # # plot(c(1:37),df[1,],type="l",main="Comparison lin vs data",xlab="Time",
# # #      ylab="value")
# # # lines(c(1:37),Pred[1,],col="red")
# # # distance<-abs(Pred-df)
# # # pred<-apply(X =abs(Pred),MARGIN = 1,FUN = calcul_norme_L2)
# # # Denominator<-apply(X = abs(df),MARGIN = 1,FUN = calcul_norme_L2)
# # # Ind_good<-which(Denominator!=0)
# # # R_squared<-(pred[Ind_good]/Denominator[Ind_good])
# # # plot(density(R_squared))
# # 
# # ###### Linear and nonlinear functions.
# # ########## PCA #############
# # matplot(t(df),type="l")
# # PCA<-FactoMineR::PCA(X = df,scale.unit = TRUE,graph = TRUE)
# # EIGEN<-PCA$svd$V
# # plot(c(1:37),EIGEN[,1])
# # plot(c(1:37),EIGEN[,2])
# # plot(c(1:37),EIGEN[,3])
# # plot(c(1:37),EIGEN[,4])
# # Quality<-PCA$ind$cos2[,1]
# # Problems_question<-as.numeric(which(Quality<quantile(Quality,0.10)))
# # MEAN_FONCTIONS<-colMeans(df[Problems_question,])
# # MEAN_opposite<-colMeans(df[-Problems_question,])
# # plot(c(1:37),y = MEAN_opposite)
# # lines(c(1:37),y=MEAN_FONCTIONS,col="red")
# # Sampled<-sample(Problems_question,size = 300)
# # matplot(t(df[Sampled,]),type="l")
# # summary(Quality)
# 
# ######### FDA ############
# # Timereg <- as.matrix((1:37)/37)
# # gaitrange  <- c(0,36)
# # Surgebasis  <- create.fourier.basis(gaitrange, nbasis=37)
# # gaitnbasis <- gaitbasis$nbasis
# # gaitcoef   <- matrix(0,gaitnbasis,dim(y)[2])
# # harmaccelLfd <- vec2Lfd(c(0, (2*pi/20)^2, 0), rangeval=gaitrange)
# # harmfdPar    <- fdPar(fd(matrix(0,37,1), Surgebasis), harmaccelLfd, 
# #                       lambda=1e5)
# # Smooth_fonctions<- smooth.basisPar(Timereg , t(df), gaitbasis,
# #                               Lfdobj=harmaccelLfd, lambda=1e-2)$fd
# # FPCA<-fda::pca.fd(Smooth_fonctions,nharm = 4,harmfdPar = harmfdPar)
# # FPCA_rotations<-varmx.pca.fd(FPCA)
# 
# # First_coeffs<-Smooth_fonctions_project$coefs[c(1:2),]
# # plot(t(First_coeffs))
# # plot.fd(Smooth_fonctions_project)
# 
# # Timereg <- as.matrix((1:36)/37)
# # gaitrange  <- c(0,35)
# # gaitbasis  <- create.fourier.basis(gaitrange, nbasis=37)
# # gaitnbasis <- gaitbasis$nbasis
# # gaitcoef   <- matrix(0,gaitnbasis,dim(y)[2])
# # harmaccelLfd <- vec2Lfd(c(0, (2*pi/35)^2, 0), rangeval=gaitrange)
# # gaitfd     <- smooth.basisPar(Timereg , t(y), gaitbasis, 
# #                               Lfdobj=harmaccelLfd, lambda=1e-2)$fd
# ##### FD Regression ######
# # Y<-gaitfd[,1]
# # X<-gaitfd[,2]
# # FDRegression<-fRegress(Y~X)
# # 
# # const  <- rep(1, dim(X$coef)[2])
# # xfdlist  <- list(const=const, Time=X)
# # 
# # beta0 <- with(Y, fd(gaitcoef, gaitbasis, fdnames))
# # beta0
# # beta1 <- with(X,  fd(gaitcoef, gaitbasis, fdnames))
# # betalist  <- list(const=fdPar(beta0), hipfd=fdPar(beta1))
# # fRegressout <- fRegress(Y, xfdlist, betalist)
# # fRegressout$betaestlist
# 
# ## Bspline basis
# 
# 
# ###################
# Nom_graph<-ifelse(Nom_v=="Surcote","S",
#                   Nom_v)
GG_violin<-ggplot(data=df$month,aes(x=factor(month_o),
                                         y=variable,
                                         fill=month_winter))+
  geom_violin()+
  geom_boxplot(alpha=0.10,col="black")+
  ylab("Surge (m)")+
  xlab("Month")+
  labs(fill="Season")+
  theme(axis.title=element_text(size=15),
        legend.text=element_text(size=10))
print(GG_violin)
# 
# 

plot(c(1:5),c(2:6),cex.lab=3,xlab="test",
     ylab="ets")
Analyse_generale<-function(Nb_sampling,period_analysis,df,l_variables,coeurs,indices_observations){

  for(i in c(1:length(l_variables))){
    if(i==1){
      nom_varibale_ill<-l_variables[i]
      obs<-parLapply(cl=coeurs,df,f_nom_variable,nom_varibale_ill)
      if(!is.null(indices_observations)){
        matrice_<-matrix(NA,nrow=length(indices_observations),ncol=37)
      }
      else{
        matrice_<-matrix(NA,nrow=length(obs),ncol=37)
      }

      for (col in c(1:37)){
        if(!is.null(indices_observations)){
          variable<-sapply(obs,FUN = function(t,serie){return(serie[t])},t=col)[indices_observations]
        }
        else{
          variable<-sapply(obs,FUN = function(t,serie){return(serie[t])},t=col)
        }
        matrice_[,col]<-variable
      }
      tirage_aleatoire<-sample(c(1:nrow(matrice_)),size = Nb_sampling,replace = FALSE)
      SS_ech<-matrice_[tirage_aleatoire,]
    }
    else{
      nom_varibale_ill<-l_variables[i]
      obs<-parLapply(cl=coeurs,df,f_nom_variable,nom_varibale_ill)
      if(!is.null(indices_observations)){
        matrice_<-matrix(NA,nrow=length(indices_observations),ncol=37)
      }
      else{
        matrice_<-matrix(NA,nrow=length(obs),ncol=37)
      }
      for (col in c(1:37)){
        if(!is.null(indices_observations)){
          variable<-sapply(obs,FUN = function(t,serie){return(serie[t])},t=col)[indices_observations]
        }
        else{
          variable<-sapply(obs,FUN = function(t,serie){return(serie[t])},t=col)
        }
        matrice_[,col]<-variable
      }
      SS_ech<-rbind.data.frame(SS_ech,matrice_[tirage_aleatoire,])
    }
  }
  colnames(SS_ech)<-c(1:ncol(SS_ech))
  re<-melt(t(SS_ech))
  colnames(re)<-c("Temps","individu","valeur")
  re$Temps<-(re$Temps-19)*(1/6)
  re$variable<-c(sapply(l_variables,function(x){return(rep(x,37*Nb_sampling))}))
  re$variable<-ifelse(re$variable=="Surcote","Surge",ifelse(re$variable=="Maree","Tide",
                                                ifelse(re$variable=="Dp.vagues.","Dp.waves",
                                                ifelse(re$variable=="Tp.vagues.","Tp.waves",
                                                ifelse(re$variable=="Dir.vent.","Dir.wind",
                                                re$variable)))))
  re$individu<-as.character((re$individu%%Nb_sampling)+1)
  GG_summary<-ggplot(re,aes(x=Temps,col=individu,group=interaction(individu),y=valeur))+
    geom_line()+
    facet_wrap(~variable,scales = "free")+
    guides(col="none")+
    xlab("Time (hour) with respect to tidal peak")+
    ylab("value")+
    theme(axis.title=element_text(size=15))
  print(GG_summary)
  Cond_one<-which(re$variable=="Surge")
  re_one<-re[Cond_one,]
  Cond_two<-which(re$variable=="Tide")
  re_two<-re[Cond_two,]
  GG_1<-ggplot(re_one,aes(x=Temps,col=individu,group=interaction(individu),y=valeur))+
    geom_line()+
    guides(col="none")+
    ylab("Surge (m)")+
    theme(axis.title=element_text(size=15))+
    theme(axis.title.x = element_blank())
  GG_2<-ggplot(re_two,aes(x=Temps,col=individu,group=interaction(individu),y=valeur))+
    geom_line()+
    guides(col="none")+
    ylab("Tide (m)")+
    theme(axis.title=element_text(size=15))+ 
    theme(axis.title.x = element_blank())
  shared_xlab <- grid::textGrob("Time (hour) with respect to tidal peak", gp = grid::gpar(fontsize = 15))
  combined<-(GG_1|GG_2) 
  resultat_graph<-patchwork::wrap_elements(full = combined) / patchwork::wrap_elements(shared_xlab) +
    plot_layout(heights = c(1, 0.05))
  print(resultat_graph)
  return(re)
}
liste_variables<-c("Surcote", "Maree","Dp.vagues.",
                   "Tp.vagues.","Dir.vent.",
                   "Hs","U")

# Import indices_donnees extrêmes -----------------------------------------
###################
Indices_chosen<-read.csv("residus_clust/inds_extremes_donnees_Surcote.csv")[,2]
TEMPS<-lapply(dataFiles,f_nom_variable,"Date")
intermed<-t(do.call(cbind.data.frame,TEMPS))
rownames(intermed)<-1:nrow(intermed)
date_compar<-intermed[,19]
colonne<-as.numeric(substr(date_compar,4,5))
vecteur_hiver<-c(1:3,9:12)
Indices_winter<-which(colonne%in%vecteur_hiver)

# apply pace --------------------------------------------------------------

saut_Pas<-3
SUb_winter<-Indices_winter[seq.int(from = 1,to = length(Indices_winter),by = saut_Pas)]
Indices_base_orig_ext<-SUb_winter[Indices_chosen]
resultat_final<-Analyse_generale(df=dataFiles,coeurs=coeurs,l_variables = liste_variables,
                                 Nb_sampling=20,period_analysis=periode,
                                 indices_observations =Indices_base_orig_ext)
resultat_final

# # Fin code ----------------------------------------------------------------
# # -----------------------------------------------------------------------
stopCluster(coeurs)
# 
