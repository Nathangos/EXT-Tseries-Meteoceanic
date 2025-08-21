require(roxygen2)

#' fonction_fichier
#'
#' @param file :str. Name of the file used.  
#'
#' @return None. Export a dataframe. 
#' @export
#'
#' @examples
fonction_fichier<-function(file){
  
  nom_colonnes<-c("Date",
                  "heure",
                  "Maree",
                  "Surcote",
                  "Hs",
                  "Tp(vagues)",
                  "Dp(vagues)",
                  "U",
                  "Dir(vent)")
  table<-file[["donnees"]]
  colnames(table)<-nom_colonnes
  nom_fichier<-substring(file[["nom"]],28,47)
  nom_dossier<-"donnees_transformees/"
  Nom_fichier_cree<-paste0(nom_dossier,nom_fichier,".csv")
  write.csv(table,file=Nom_fichier_cree)
}

fonc_skip<-function(file){
  return(list("donnees"=read.table(file,skip=9),"nom"=file))
}

# Pareto definitions ------------------------------------------------------

#' get.tail.index.estimate
#'
#' @param vec : vector. Non-extreme and extreme observations.
#' @param p : float. Weight put on the bulk of the density.
#' @param sigma.fixed : float. Estimated scale parameter. 
#'
#' @return estimated ML shape parameter. 
#' @export
#'
#' @examples
get.tail.index.estimate <-function(vec, p, sigma.fixed){
  
  thresh=quantile(vec,p)-10^{-10}
  is.exc=which(vec-thresh>0)
  fun2opt=function(par){
    -sum(log(dgp_craft(x=vec[is.exc]-thresh,sigma=sigma.fixed,par)))
  }
  optim(par=0,fun2opt,method="Brent",lower=-1,upper=5)$par
}

#' dgp_craft
#'
#' @param x : float. Extreme or non-extreme observations. 
#' @param sigma : float. Estimated scale parameter. 
#' @param xi : float. Estimated shape parameter. 
#'
#' @return GPD density obtained for the x value. 
#' @export
#'
#' @examples
dgp_craft<-function(x,sigma,xi){
  
  if(abs(xi)<10^{-4}){ # avoid numerical instabilities
    exp(-x/sigma)/sigma
  }else{
    result = ifelse(1+xi*x/sigma<=0,0,sigma^{-1}*(1+xi*x/sigma)^{-1/xi-1})
    ifelse(is.na(result) | (result == 0), .Machine$double.xmin, result)
  }
}

#' qgp_craft
#'
#' @param x : float. Desired level of the quantile. 
#' @param sigma : float. Estimated scale parameter. 
#' @param xi : float. Shape parameter. 
#'
#' @return Quantile of a level of the chosen level. 
#' @export
#'
#' @examples
qgp_craft<-function(x,sigma,xi){
  
  if(abs(xi)<10^{-4}){ # avoid numerical instabilities
    -sigma * log(1-x)
  }else{
    sigma/xi *((1-x)^{-xi}-1)
  }
}

#' pgp_craft
#'
#' @param x: float. Extreme or non-extreme value 
#' @param sigma: float. Estimated scale parameter.  
#' @param xi: float. Estimated shape parameter. 
#'
#' @return Probability put on the input value. 
#' @export
#'
#' @examples
pgp_craft<-function(x,sigma,xi){
  
  if(abs(xi)<10^{-4}){ # avoid numerical instabilities
    1-exp(-x/sigma)
  }else{
    1- ifelse(1+xi*x/sigma<=0, 0, (1+xi*x/sigma)^{-1/xi})
  }
}

# Qualite transformation --------------------------------------------------
#' test_Uniforme
#'
#' @param collection 
#'
#' @return p_valeur du test d'adequation a une uniforme(-1,0).
#' @export
#'
#' @examples
test_Uniforme<-function(collection_dimension){
  KSTESTunif<-ks.test(collection_dimension,"punif",min=-1,max=0)
  return(KSTESTunif$p.value)
}
##### Calculer la norme de chaque série temporelle. #####

#' calcul_norm_L2
#'
#' @param series: vector(series). Time series. 
#'
#' @return float. L2 norm of the series (positivity considered).
#' @export
#'
#' @examples
calcul_norm_L2<-function(series){
  L<-length(series)
  pas<-(1/(L-1))

  # Positivity indicator -----------------------------------------------
  Sigma<-sum(subset(series,series>0)**2)
  norme_carree<-Sigma*pas
  return(norme_carree**(1/2))
}

d_cste_norme_L2<-function(series){
  L<-length(series)
  pas<-1/L
  Sigma_inv<-pas^(1/2)*(sum(subset(series,series>0)**2))^(-1/2)
  return(Sigma_inv)
}

fdiff<-function(series){
  return(diff(series))
}
#' fonction_calcul_norme
#'
#' @param df 
#' @param nom_colonne : column used to calculate the norm. 
#'
#' @return float.
#' @export
#'
#' @examples
fonction_calcul_norme<-function(df,nom_colonne){
  
  series<-df[,nom_colonne]
  resultat_calcul_norme<-calcul_norm_L2(series)
  return(resultat_calcul_norme)
}

#' f_nom_par_variable
#'
#' @param file : file name. 
#' @param nom_var : variable name to return
#'
#' @return vector. Column of the variable. 
#' @export
#'
#' @examples
f_nom_variable<-function(file,nom_var){
  
  return(file[,nom_var])
}
#' fonction_renvoi_angle
#'
#' @param df 
#' @param nom_colonne 
#'
#' @return the angular component of a series. 
#' @export
#'
#' @examples
fonction_renvoi_angle<-function(donnees){
  
  norme<-calcul_norm_L2(donnees)
  if (norme>0){
    angle<-donnees/norme
    return(angle)
  }
  else{
    return(donnees)
  }
  
}

####### Obtenir les résidus d'un modèle ARIMA. ####
fonction_resid_arima<-function(series,p,d,q){
  modele<-arima(x = series,order=c(p,d,q))
  return(modele$residuals)
}
#### Observer l'évolution pour chaque coordonnée au cours des années.###

fonction_reg_par_coord<-function(colonne_coord){
  nb_obs<-length(colonne_coord)
  x_reg<-c(1:nb_obs)
  reg_lin<-lm(colonne_coord~x_reg)
  return(reg_lin)
}

#### Hill. ####
#' function_estimator_hill
#'
#' @param k : rank in the descreasing order of the statistic. 
#' @param series_ ; series of observation. 
#'
#' @return Float. Shape parameter of the distribution of the maximum/excesses. 
#' @export
#'
#' @examples
function_estimator_hill<-function(k,series_){
  
  ord_series_norm<-order(series_,decreasing = TRUE)
  series_norm_sort<-series_[ord_series_norm]
  dernier<-(k-1)
  U_ndk<-series_norm_sort[k]
  ss_series_norm_sort<-series_norm_sort[c(1:dernier)]
  gamma_hill<-mean(log(ss_series_norm_sort/U_ndk))
  return(gamma_hill)
}
##### Calculate the moments estimator. #####
#' function_estimateur_moment
#'
#' @param k 
#' @param series_ 
#'
#' @return Float. Moment estimator of the shape parameter. 
#' @export
#'
#' @examples
function_estimator_moment<-function(k,series_){
  
  ord_series_norm<-order(series_,decreasing = TRUE)
  series_norm_sort<-series_[ord_series_norm]
  dernier<-(k-1)
  U_ndk<-series_norm_sort[k]
  ss_series_norm_sort<-series_norm_sort[c(1:dernier)]
  gamma_hill<-mean(log(ss_series_norm_sort/U_ndk))
  M2<-mean((log(ss_series_norm_sort/U_ndk))^2)
  M1<-max(gamma_hill,0)
  gamma_inf<-1-(1/2)*(1-((M1)^2/M2))^(-1)
  estimator_moment<-M1+gamma_inf
  return(estimator_moment)
}
POT_fonction_estimateur_moment<-function(k,series){
  seuil<-sort(series,decreasing = TRUE)[k]
  modele<-POT::fitgpd(data = series,threshold = seuil)
  return(as.numeric(modele$param[2]))
  
}

# Estimateur_Likelihood ---------------------------------------------------
#' fonction_ML_extRemes
#'
#' @param k 
#' @param data_d
#' @param typeML
#'
#' @return extRemes package estimator of the shape parameter with its boundaries. 
#' @export
#'
#' @examples
function_ML_extRemes<-function(k,data_d,typeML="GP",NB_years=NULL){
  
  ## parametric method
  data_ord<-sort(data_d,decreasing = TRUE)
  th<-data_ord[k]
  L<-length(data_d)/NB_years
  if(is.null(NB_years)){
    theta_k<-extRemes::fevd(data_d,threshold = th,type=typeML)
  }
  else{
    theta_k<-extRemes::fevd(data_d,threshold = th,type=typeML,
                            time.units = paste0(L,"/year"))
  }
  
  tryCatch(expr={
    r<-distillery::ci(theta_k,type = "parameter")
    # formula is
    # std_<-sqrt(inverse_hessian)[2,2]
    # ga-qnorm(0.975)*std_
    # ga+qnorm(0.975)*std_
    q_minus<-r[2,1]
    estimator<-r[2,2]
    q_plus<-r[2,3]
    liste_r<-list("q_minus"=q_minus,"q_plus"=q_plus,"estimator"=estimator)
    return(liste_r)
    },
    error=function(e){
      estimator<-theta_k$results$par[[2]]
      liste_r<-list("q_minus"=NA,"q_plus"=NA,"estimator"=estimator)
      return(liste_r)}
  )
}
#' fonction_MLplot_resume
#'
#' @param resultatML :result from fonction_ML_extRemes
#' @param vecteur_k : vector(int). Vector of the order used. 
#' @param nom_variable : str. Name of the variable analyzed. 
#'
#' @return Ggplot.
#' @export
#'
#' @examples
fonction_MLplot_resume<-function(resultatML,vecteur_k,nom_variable,lims_Y=c(-2,2)){
  
  intermed<-t(resultatML)
  mink<-min(vecteur_k)
  maxk<-max(vecteur_k)
  bilan<-as.data.frame(apply(intermed,MARGIN = 2,FUN=unlist))
  nom_plot<-paste("Evolution du gamma pour",nom_variable)
  GA_plot<-ggplot(data=bilan,aes(x=vecteur_k,y=estimator,col="estimateur"))+
    geom_line()+
    ggtitle(nom_plot)
  if (all(!is.na(bilan))==TRUE){
    GA_plot<-GA_plot+
      geom_line(aes(x=vecteur_k,y = q_minus,col="borne inferieure"))+
      geom_line(aes(x=vecteur_k,y = q_plus,col="borne superieure"))+
      ylim(lims_Y)
  }
  GA_plot<-GA_plot+
    labs(caption=paste("k going from",mink,"to",maxk),colour="Legend",x="Valeur de k")

  return(GA_plot) 
}

# Ajouter points/reconstruction -------------------------------------------
f_approx_interpol<-function(series)
{
  nb<-length(series)
  vecteur_temps<-c(1:nb)*1/nb
  Interpol_resultat<-approx(vecteur_temps,series)
  return(Interpol_resultat$y)
}

#### Calculer les scores ###

#' fonction_score_distrib
#'
#' @param rang_fonction : the rank of the eigen function in the contribution to the variance. 
#' @param matrice_fct_propre : the set of eigen functions. 
#' @param Mat_obs_centrees : the set of centered observations. 
#'
#' @return vector(float). Coordinates of the observations on the dimension.
#' @export
#'
#' @examples
fonction_score_distrib<-function(rang_fonction,matrice_fct_propre,Mat_obs_centrees){
  
  Fonction_prop<-matrice_fct_propre[,rang_fonction]
  Fonction_prop
  Scores_trouves<-as.matrix(Mat_obs_centrees)%*% unlist(Fonction_prop)
  return(Scores_trouves)
}

#### Simulation_nouveau_score ###

#' fonction_simul_nouveau_score
#'
#' @param variable_uniforme : result from a runif.
#' @param fonction_repart : result of the ecdf function.
#' @param ech : vector(float). Sample of coordinates. 
#'
#' @return Float. The coordinate sampled. 
#' @export
#'
#' @examples
fonction_simul_nouveau_score<-function(variable_uniforme,fonction_repart,ech){
  
  valeurs_fonct_points<-fonction_repart(ech)
  distance<-abs(valeurs_fonct_points-variable_uniforme)
  return(max(ech[which(distance==min(distance))]))
}
#' fonction_simulation_gpareto
#'
#' @param loc : location parameter. 
#' @param scale : scale parameter (>0)
#' @param shape : shape parameter.
#'
#' @return Float. Simulate a realization of a GP. 
#' @export
#'
#' @examples
fonction_simulation_gpareto<-function(loc,scale,shape){
  
  fait<-FALSE
  # contrainte sur (1+ gamma/scale x) à prendre en compte.
  while (fait==FALSE){
    rapport <-(shape/scale)
    u<-runif(1,min = 0,max=1)
    x<-((1-u)^((-1)*shape)-1)*(scale/shape)+loc
    if (1+(x-loc)*rapport>0){
      fait=TRUE
    }
  }
  return(x)
}

#' f_marginales_uniforme
#'
#' @param variable : float. The sample of observations used. 
#' @param p_u : float. The threshold used to distinguish extreme values. 
#' @param n.dens : int. The number of points used to define the kernel density estimator. 
#'
#' @return Vector(float). Sample of transformed observations.
#' @export
#'
#' @examples
f_marginales_uniforme<-function(variable,p_u,n.dens){
  
  kernel_dens<-density(x=variable,n=n.dens)
  seuil<-as.numeric(quantile(variable,1-p_u))
  extremes<-subset(variable,variable>seuil)
  indice<-which.max(seuil-kernel_dens$x<0)-1
  valeur_fonc_u<-as.numeric(kernel_dens$y[indice])
  scale_fonc<-(p_u/valeur_fonc_u)
  #indice pour k. 
  k<-which.max(variable-seuil<0)-1
  theta<-extRemes::fevd(variable,threshold = seuil,type = "GP")
  gamma<-theta$results$par[[2]]
  # indices_moins
  ecart<-diff(kernel_dens$x)[1]
  # Riemann
  F_dens<-cumsum(kernel_dens$y)*ecart
  non_ext<-subset(variable,variable<=seuil)
  intervals<-c(-Inf,kernel_dens$x+ecart/2)
  interval_vecteur<-findInterval(non_ext,vec = intervals)
  # nouvelles valeurs
  n<-length(variable)
  vecteur_unif<-rep(NA,length(variable))
  vecteur_unif[variable<=seuil]<-F_dens[interval_vecteur]
  vecteur_unif[variable>seuil]<-(1-p_u)+(p_u)*evd::pgpd(extremes,loc=seuil,scale=scale_fonc,shape=gamma)
  vecteur_unif_decalage<-vecteur_unif-1
  return(as.vector(vecteur_unif_decalage))
}

#' fonc_norm_inv
#'
#' @param variable : vector(float).
#'
#' @return Float. 
#' @export
#'
#' @examples
fonc_norm_inv<-function(variable){
  vecteur_variable_transformee<-as.numeric(-(variable)^(-1))
  valeur_norme<-calcul_norm_L2(vecteur_variable_transformee)
  valeur_norme_geo<-(-1/valeur_norme)
  return(valeur_norme_geo)
}
#' f_NWatson
#'
#' @param vecteur_x : vector(float). Time coordinates ("training" sample).
#' @param x : float. The new input. 
#' @param vecteur_y : vector(float). Output values ("training" sample)
#' @param h : float. Bandwidth. 
#'
#' @return Float. Value of the Nadaraya-Watson estimator. 
#' @export
#'
#' @examples
f_NWatson<-function(vecteur_x,x,vecteur_y,h){
  
  numerateur<-sum(vecteur_y*dnorm((vecteur_x-x/h)))
  denominateur<-sum(dnorm((vecteur_x-x)/h))
  return(numerateur/denominateur)
}
#' f_lissage
#'
#' @param donnees 
#' @param M : int. Number of replicas. 
#' @param h : float. The value of the bandwidth. 
#'
#' @return Vector(float). Results from N-Watson estimator. 
#' @export
#'
#' @examples
f_lissage<-function(donnees,M,h){
  
  L<-length(donnees)
  vecteur_x<-seq.int(from = 0,to=1,length.out = L)
  echantillon_plus_fin<-seq.int(from=0,to=1,length.out = M)
  nouvelles_donnees<-sapply(echantillon_plus_fin,FUN = f_NWatson,vecteur_x=vecteur_x,vecteur_y=donnees,h=h)
  return(nouvelles_donnees)
}
fonc_Median_inv<-function(x){
  inverse<--1/x
  R<--1/median(inverse)
  return(R)
}
f_quantile<-function(donnees,l_tau){
  P<-ncol(donnees)
  liste_t<-c(1:P)
  vecteur_p_u<-rep(NA,P)
  if (length(l_tau)>1){
    for (t in liste_t){
      variable_t<-unlist(donnees[,t])
      tau<-l_tau[t]
      quantile<-mean(variable_t<=l_tau[t])
      vecteur_p_u[t]<--quantile+1
    }
  }
  else{
    tau<-l_tau
    for (t in liste_t){
      variable_t<-unlist(donnees[,t])
      quantile<-mean(variable_t<=tau)
      vecteur_p_u[t]<--quantile+1
    }
  }
  return(vecteur_p_u)
}
fonction_quantile_tau<-function(coeurs,nom_variable,df,opt_diff,liste_t,l_tau){
  print(paste("On utilise comme u moyen,",mean(l_tau),"avec une differenciation=",opt_diff))
  obs<-parLapply(cl=coeurs,df,f_nom_variable,nom_variable)
  obs<-do.call(rbind.data.frame,obs)
  if(opt_diff==TRUE){
    obs<-as.data.frame(t(diff(t(as.matrix(obs)))))
    colnames(obs)<-c(1:length(colnames(obs)))
  }
  colnames(obs)<-c(1:length(colnames(obs)))
  vecteur_pu<-f_quantile(donnees=obs,l_tau=l_tau)
  return(vecteur_pu)
}


# Fonction_qualite ajustement ---------------------------------------------
APPEL_KS_AD<-function(indice,donnees,l_seuil){
  seuil<-l_seuil[indice]
  echantillon<-donnees[,indice]
  return(fonction_KSAD_GP(echantillon,seuil))
}
fonction_KSAD_GP<-function(echantillon,seuil){
  Estimateurs<-extRemes::fevd(echantillon,threshold = seuil,type="GP",method = "MLE")$result$par
  Scale<-Estimateurs[[1]]
  Shape<-Estimateurs[[2]]
  excedents<-subset(echantillon,echantillon>seuil)
  KS<-ks.test(excedents,extRemes::"pevd",threshold=seuil,scale=Scale,shape=Shape,type="GP")$p.value
  AD<-goftest::ad.test(excedents,null=extRemes::"pevd",threshold=seuil,scale=Scale,shape=Shape,type="GP",estimated = FALSE)$p.value
  # prendre en compte le fait qu'on utilise des estimateurs et pas les 
  # vraies valeurs. 
  AD_estim<-goftest::ad.test(excedents,null=extRemes::"pevd",threshold=seuil,scale=Scale,shape=Shape,type="GP",estimated = TRUE)$p.value
  return(list("KSmirnov"=KS,"ADarling"=AD,"AD_BRAUN"=AD_estim))
}
#' Generation_Pareto_std
#'
#' @param u_pareto : float. Realisation of the uniform law
#'
#' @return Realisation of a standard Pareto law. 
#' @export
#'
#' @examples
Generation_Pareto_std<-function(u_pareto){
  
  return(1/(1-u_pareto))
}
fnct_cov_par_temps<-function(t,s,f_variogram,f_sigma){
  cov_ts<-(f_sigma(h=t)+f_sigma(h=s))-f_variogram(t=t,s=s)
  return(cov_ts)
}


Generation_Sigma<-function(f_variogram,f_sigma,Vecteur_temps){
  L<-length(Vecteur_temps)
  vecteur_zero<-rep(0,L*L)
  Matrice_cov_moitie<-matrix(vecteur_zero,byrow = TRUE,nrow=L)
  for(j in 1:L){
    vecteur_autres_temps<-Vecteur_temps[c(1:j)]
    Ligne<-sapply(vecteur_autres_temps,FUN=fnct_cov_par_temps,f_variogram=f_variogram,f_sigma=f_sigma,s=Vecteur_temps[j])
    Matrice_cov_moitie[j,c(1:j)]<-Ligne    
  }
  Matrice_cov_reste<-t(Matrice_cov_moitie)
  Matrice_cov<-Matrice_cov_moitie+Matrice_cov_reste
  diag(Matrice_cov)<-diag(Matrice_cov)/2
  return(Matrice_cov)
}
Trajectoire_gp<-function(f_variogram,f_sigma,Vecteur_temps,nb,Parametres){
  coord <- Vecteur_temps
  Sigma<-Generation_Sigma(f_variogram = f_variogram,f_sigma=f_sigma,Vecteur_temps = Vecteur_temps)
  Sigma_racine_caree<-chol(Sigma)
  Realisations_GP<-mvtnorm::rmvnorm(n=nb,sigma = Sigma,method="chol")
  return(Realisations_GP)
}
Log_norm_transf<-function(une_traject_GP,f_sigma,param_alpha,Vecteur_temps){
  liste_sigma<-sapply(Vecteur_temps,f_sigma)
  log_norm_process<-exp((1/param_alpha)*(une_traject_GP-(liste_sigma)))
  return(log_norm_process)
}

# Produce conditional simul of BR--simul of standard residuals -----------------------------------------
#' fonction_simul_BR_accept_reject
#'
#' @param NB_BR : desired number of simulations
#' @param modele_BR : output of fitmaxstab for the BR model
#' @param seuil_Frech : threshold used to identify extreme residuals 
#'
#' @return Df of standard residuals.
#' @export
#'
#' @examples
fonction_simul_BR_accept_reject<-function(NB_BR,modele_BR,seuil_Frech){
  J<-0
  Sorties<-matrix(NA,nrow = NB_BR,ncol = 37)
  deb<-1
  while(J<NB_BR){
    Simul<-SpatialExtremes::rmaxstab(NB_BR-J,Locations,cov.mod="brown",
                                     range=modele_BR$param[1], 
                                     smooth=modele_BR$param[2])
    L2<-apply(X = Simul,MARGIN = 1,FUN = calcul_norm_L2)
    inds_taken<-which(L2<seuil_Frech)
    end<-length(inds_taken)+deb-1
    indices_base<-c(deb:end)
    Sorties[indices_base,]<-Simul[inds_taken,]
    deb<-end+1
    J<-end
  }
  return(Sorties)
}

Trajectoire_log_norm<-function(f_variogram,f_sigma,Vecteur_temps,nb,alpha){
  GP<-Trajectoire_gp(f_variogram=f_variogram,Vecteur_temps=Vecteur_temps,nb=nb,f_sigma=f_sigma)
  LOG_Norm_Ensemble<-apply(GP,MARGIN=1,FUN =Log_norm_transf,f_sigma=f_sigma,param_alpha=alpha,Vecteur_temps=Vecteur_temps)
  return(LOG_Norm_Ensemble)
}
#' Procedure_MHastings
#'
#' @param echantillons_log_norm : Dataframe. Trajectoires d'un processus log-normal.
#' @param Longueur_echantillon 
#'
#' @return list. M Processus de forme W/l(W). 
#' @export
#'
#' @examples
Procedure_MHastings<-function(echantillons_log_norm,Longueur_echantillon){
  
  L<-nrow(echantillons_log_norm)
  liste_sigma_l<-list()
  premiere_traject<-echantillons_log_norm[1,]
  premiere_norme<-calcul_norm_L2(premiere_traject)
  liste_realisation_MH<-list()
  liste_realisation_MH[[1]]<-premiere_traject
  liste_lforme<-list()
  liste_lforme[[1]]<-(premiere_traject/premiere_norme)
  Q<-c()
  Moy_Moy<-c()
  Moy_Max<-c()
  for(i in (2:L)){
    passe<-(i-1)
    #MH a verifier. L'expression est différente de ce qu'on voit d'habitude. 
    rapport<-(calcul_norm_L2(echantillons_log_norm[i,])/calcul_norm_L2(series = liste_realisation_MH[[passe]]))
    p_n<-min(rapport,1)
    if(p_n==1){
      U1<-1
    }
    else{
      u<-runif(n = 1)
      U1<-as.numeric(u<=p_n)
    }
    if(U1==0){
      Q<-c(Q,0)
      trajectoire_conservee<-liste_realisation_MH[[passe]]
    }
    else{
      Q<-c(Q,1)
      trajectoire_conservee<-echantillons_log_norm[i,]
    }
    Norme<-calcul_norm_L2(trajectoire_conservee)
    liste_realisation_MH[[i]]<-trajectoire_conservee
    candidat_theta<-(trajectoire_conservee/Norme)
    Moy<-mean(candidat_theta)
    Max_conv<-max(candidat_theta)
    Moy_Max<-c(Moy_Max,Max_conv)
    Moy_Moy<-c(Moy_Moy,Moy)
    liste_lforme[[i]]<-candidat_theta
  }
  # graphique convergence ---------------------------------------------------
  plot(cummean(Q),type="l",xlab="Iteration",ylab="Moyenne cumulee",main=paste0("Evolution de la probabilité d'acceptation pour ",L," simulations"))
  plot(cummean(Moy_Moy),type="l",xlab="Iteration",ylab="Moyenne cumulee",main=paste0("Evolution du niveau moyen pour ",L," simulations" ))
  plot(cummean(Moy_Max),type="l",xlab="Iteration",ylab="Moyenne cumulee",main=paste0("Evolution du maximum pour ",L," simulations"))
  indices_pris<-(L-Longueur_echantillon)+1
  indice_fin<-L
  return(liste_lforme[c(indices_pris:indice_fin)])
}

variogram_alpha_lambda<-function(ALPHA,lambda,t,s){
  h<-(t-s)
  rapport<-(abs(h)/lambda)
  semi_vario<-rapport^(ALPHA)
  return(semi_vario)
}
sigma_alpha_lambda<-function(ALPHA,lambda,t){
  rapport<-(t/lambda)
  semi_vario<-rapport^(ALPHA)
  return(semi_vario)
}

estim_Pearson_tidal_cycle<-function(Matrice_couples,individus_select,type_corr){
  Pearson_empirique<-sapply(Matrice_couples,FUN=Pearson_t_s,
                            functional_inds=individus_select,
                            type_corr=type_corr)
  return(Pearson_empirique)
}
Pearson_t_s<-function(functional_inds,couples_indices,type_corr){
  t<-couples_indices[1]
  s<-couples_indices[2]
  return(cor(x = functional_inds[,t],y=functional_inds[,s],
             method = type_corr))
}

inv_extremogram_empirique<-function(Matrice_couples,individus_select,Tau){
  PIST_empirique<-sapply(Matrice_couples,FUN = pi_s_t,functional_inds=individus_select,vecteur_u=Tau)
  return(PIST_empirique)
}
pi_s_t<-function(functional_inds,vector_u,couples_indices){
  t<-couples_indices[1]
  s<-couples_indices[2]
  if(is.null(nrow(functional_inds))==FALSE){
    numerateur<-sum(as.numeric((functional_inds[,s]>vector_u[s])&(functional_inds[,t]>vector_u[t])))
    denominateur<-sum(as.numeric(functional_inds[,s]>vector_u[s]))
  }
  else{
    numerateur<-sum(as.numeric((functional_inds[s]>vector_u[s])&(functional_inds[t]>vector_u[t])))
    denominateur<-sum(as.numeric(functional_inds[s]>vector_u[s]))
  }
  return(numerateur/denominateur)
}

inv_pi_s_t<-function(functional_inds,vecteur_u,couples_indices){
  t<-couples_indices[1]
  s<-couples_indices[2]
  numerateur<-sum(as.numeric((functional_inds[,s]<=vecteur_u[s])&(functional_inds[,t]<=vecteur_u[t])))
  denominateur<-sum(as.numeric(functional_inds[,s]<=vecteur_u[s]))
  return(numerateur/denominateur)
}
#' function_reconstitution_trajectory_std
#'
#' @param Vector_coords : dataframe. Pca coordinates obtained.
#' @param Base_functions_p : dataframe. Eigen vectors obtained. 
#' @param NB_dim : int. Chosen number of PCA vectors. 
#' @param mu_t : vector(float). Mean function. 
#' @param var_t : vector(float). Standard value function. 

#'
#' @return Angular components of the observations. 
#' @export
#'
#' @examples
function_reconstitution_trajectory_std<-function(Vector_coords,Base_functions_p,NB_dim,mu_t,sd_t){
  L<-nrow(Vector_coords)
  Matrix_score<-diag(Vector_coords[,1],ncol = L)
  feigen1<-Base_functions_p[,1]
  Matrice_repmfonction<-matrix(rep(feigen1,L),ncol=L)
  Simulations<-t(Matrice_repmfonction%*%Matrix_score)
  Mat_mean<-matrix(rep(mu_t,L),byrow = TRUE,nrow=L)
  if(NB_dim==1){
    Shape_ACP<-Simulations%*%diag(sd_t)+Mat_mean
  }
  else{
    for(j in c(2:NB_dim)){
      MS<-diag(Vector_coords[,j] ,ncol=L)
      Matrice_f_p2<-matrix(rep(Base_functions_p[,j],L),ncol=L)
      Simulations<-Simulations+t(Matrice_f_p2%*%MS)
    }
    Shape_ACP<-Simulations%*%diag(sd_t)+Mat_mean
  }
  
  return(Shape_ACP)
}

#' Title
#'
#' @param Vector_coords 
#' @param Base_functions_p 
#' @param NB_dim 
#' @param mu_t 
#'
#' @return Composantes angulaires des observations. 
#' @export
#'
#' @examples
fonction_reconstitution_trajectoire<-function(Vecteur_coords,Base_fonctions_p,NB_dim,mu_t){
  L<-nrow(Vecteur_coords)
  Matrice_score<-diag(Vecteur_coords[,1],ncol = L)
  fpropre1<-Base_fonctions_p[,1]
  Matrice_repmfonction<-matrix(rep(fpropre1,L),ncol=L)
  Simulations<-t(Matrice_repmfonction%*%Matrice_score)
  Mat_moyenne<-matrix(rep(mu_t,L),byrow = TRUE,nrow=L)
  if(NB_dim==1){
    FORME_ACP<-Simulations+Mat_moyenne
  }
  else{
    for(j in c(2:NB_dim)){
      MS<-diag(Vecteur_coords[,j] ,ncol=L)
      Matrice_f_p2<-matrix(rep(Base_fonctions_p[,j],L),ncol=L)
      Simulations<-Simulations+t(Matrice_f_p2%*%MS)
    }
    # Deuxieme score ----------------------------------------------------------
    FORME_ACP<-Simulations+Mat_moyenne
  }
  
  return(FORME_ACP)
}

#' function_Structure_Matrix
#'
#' @param NB_dim : int. Number of coordinates
#'
#' @return Generate the matrix used in Rvine for the structure
#' @export
#'
#' @examples
function_Structure_Matrix<-function(NB_dim){
  
  Matrix_C <- c()
  for(j in c(1:NB_dim)){
    lim<-NB_dim-1
    vector<-rep(0,NB_dim)
    if(j<lim){
      vector[1:j]<-rep(j,j)
    }
    if(j==1){
      vector[1]<-NB_dim
    }
    if(j==lim){
      vector[1:j]<-rep(j,j)
      vector[1]<-1
    }
    if(j==NB_dim){
      vector[1]<-lim
      vector[2:NB_dim]<-rep(1,lim)
    }
    Matrix_C<-c(Matrix_C,vector)
  }
  
  Matrix_C <- matrix(Matrix_C, NB_dim, NB_dim,byrow = TRUE)
  return(Matrix_C)
}


FNCT_percentile_ppourcent<-function(series,debut,pas){
  vecteur_quantile<-c()
  for(j in c(1:3)){
    vecteur_quantile<-c(vecteur_quantile,as.numeric(quantile(series,debut+j*pas)))
  }
  vecteur_quantile<-c(vecteur_quantile,as.numeric(quantile(series,0.99)))
  ecart_type<-sd(series)
  mu_t<-mean(series)
  return(c(mu_t,ecart_type,vecteur_quantile))
}
###############
#' calcul S_theta
#'
#' @param ALPHA (float): parametre
#' @param LAMBDA (float) : parametre
#' @param L (int) : pas des observations
#' @param individus_select (df) : individus selectiones
#' @param Tau (vecteur) : vecteur des seuils GP. 
#'
#' @return float. Somme des erreurs au carre entre théorie et réalité.
#' @export
#'
#' @examples
fonction_Stheta<-function(theta,Matrice_couples,L,PIST_empirique,fonction_theorique){
  vecteur_temps<-(1:L)/L
  #Theorique
  if(theta[2]>2){return(1e50)}
  f_variog<-function(x,y){
    h<-abs(x-y)
    return(fonction_theorique(h,theta))
  }
  PIST_theorique<-sapply(Matrice_couples,FUN=pi_s_t_theorique,vecteur_temps=vecteur_temps,f_variogram=f_variog)
  S_theta<-mean((PIST_theorique-PIST_empirique)^(2))
  return(S_theta)
}
pi_s_t_theorique<-function(couples_indices,vecteur_temps,f_variogram){
  t<-vecteur_temps[couples_indices[1]]
  s<-vecteur_temps[couples_indices[2]]
  valeur_aleatoire<-(f_variogram(t,s)/2)^(1/2)
  pist_theorique<-2*(1-pnorm(valeur_aleatoire))
  return(pist_theorique)
}

fonction_creation_matrice<-function(M){
  vecteur<-1:M
  fcnt_each<-function(coord_une,vecteur_coord){
    return(sapply(vecteur_coord,function(c1,c2){return(c(c1,c2))},c1=coord_une))
  }
  return(lapply(vecteur,fcnt_each,vecteur_coord=vecteur))
}
fnct_select_colonne<-function(df,numero){
  return(df[numero,])
}
fonction_trajectoire_positive<-function(individu){
  resultat_test<-all(individu>0)
  return(resultat_test)
}

# fonctions_MV_clust ------------------------------------------------------

fnct_simul_one<-function(vecteur_sigma,vecteur_mu,prob_classe){
  u<-rbinom(n = 1,size = 1,prob=prob_classe[1])
  if(u==1){
    realisation<-rnorm(1,mean=vecteur_mu[1],sd = vecteur_sigma[1])
  }
  else{
    realisation<-rnorm(1,mean=vecteur_mu[2],sd = vecteur_sigma[2])
  }
  return(realisation)
}
fnct_simul_MV_MCLUST<-function(objetMCLUST,M){
  PARAMS_MC<-objetMCLUST$parameters
  vecteur_mu<-PARAMS_MC$mean
  probabilite_classes<-PARAMS_MC$pro
  vecteur_sigma<-sapply(PARAMS_MC$variance$sigmasq,sqrt)
  realisations<-replicate(M,expr = fnct_simul_one(vecteur_sigma = vecteur_sigma,vecteur_mu = vecteur_mu,prob_classe = probabilite_classes))
  return(realisations)
}
Qpareto<-function(x){
  numerateur<-(1-x)^(-1)-1
  denom<-1
  y<-((numerateur/1)*1)+1
  return(y)
}
Qfrech<-function(x){
  return(-log(x)^(-1))
}
graphique_qlog<-function(series_base,series_simulations,nom_variable,debut_prop,origine_simul){
  probabilites<-seq.int(from = debut_prop,to = 0.99,by = 0.01)
  obj_qqplot<-extRemes::qqplot(series_base,series_base,make.plot=FALSE)
  debut<-obj_qqplot$qdata$x
  upper<-obj_qqplot$qdata$upper
  lower<-obj_qqplot$qdata$lower
  indices<-which((is.na(lower)==FALSE)&(is.na(upper)==FALSE))
  # probabilites autre base. 
  
  # f_emp<-ecdf(x =series_simulations)
  # P<-sort(f_emp(series_simulations))
  # ind_tronc<-which(P>=debut_prop)
  estimateur<-as.numeric(quantile(debut,probabilites))
  up_quantile<-as.numeric(quantile(upper[indices],probabilites))
  low_quantile<-as.numeric(quantile(lower[indices],probabilites))
  Min_logy<-min(c(min(low_quantile),min(up_quantile),min(estimateur),min(series_simulations)))
  Max_logy<-max(c(max(low_quantile),max(up_quantile),max(estimateur),max(series_simulations)))
  series_simul2<-as.numeric(quantile(series_simulations,probabilites))
  # plot(probabilites,estimateur,log="xy",type="o",ylab=paste0(nom_variable," (échelle log)"),ylim=c(Min_logy,Max_logy),main=paste0("QQplot de la ",nom_variable," (log-log)"),xlab="probabilité de non-dépassement",col=1)
  # points(probs_trouvees,series_simul2,col=2,type="l")
  # lines(probabilites,up_quantile,lty=2,col=3)
  # lines(probabilites,low_quantile,lty=2,col=3)
  # legend("topleft",legend = c("Observations","Bandes de confiance (KS)",paste0("Simulations (",origine_simul,")")),col = c(1,2,3),lty=c(1,1,2),title = "Variable",cex=0.8)
  resume<-list(probs=probabilites,"Percentile_95_KS"=up_quantile,"Percentile_5_KS"=low_quantile,"Simulations"=series_simul2,"estimateur_données"=estimateur,"origine"=origine_simul)
  return(resume)
}

#' Title
#'
#' @param matrice_1 : premiere matrice
#' @param matrice_2 : seconde matrice. 
#' @param nom_variable : str. Nom de la fonctionnelle étudiée. 
#'
#' @return Matrice des p valeurs de Mann-Whitney et de KS. 
#' @export
#'
#' @examples
f_p_valeur_KS_MW<-function(matrice_1,matrice_2,nom_variable, nom_methode){
  
  vecteur_colonnes<-c(1:ncol(matrice_1))
  resultat<-sapply(X = vecteur_colonnes,function(t,Mat,Mat2){
    mwt<-wilcox.test(Mat[,t],Mat2[,t])
    kst<-ks.test(Mat[,t],Mat2[,t])
    return(list("KS"=kst$p.value,"MW"=mwt$p.value))},Mat=matrice_1,Mat2=matrice_2)
  Df<-as.data.frame(t(resultat))
  Df<-melt(sapply(Df,unlist))
  colnames(Df)<-c("Temps","Type","p.valeur")
  GGmarg<-ggplot2::ggplot(data=Df,aes(x = Temps,y = p.valeur,col=Type,group=interaction(Type)))+
    geom_point()+
    geom_line()+
    ylim(0,1)+
    geom_hline(yintercept = 0.05,show.legend = TRUE)+
    ggtitle(paste0("Test des marginales pour la variable (",nom_variable,") avec ",nom_methode))
  print(GGmarg)
}

#' Title
#'
#' @param base_RV : dataframe respectant les variations régulières. 
#' @param nb_scores : nombre de scores (PCA) à étudier. 
#' @param stand_ : standardisation utilisée ou non
#'
#' @return
#' @export
#'
#' @examples
Analyse_extreme_proj<-function(base_RV,nb_scores,stand_,nom_variable,L,NB_year,M1,M2){
  ACP<-FactoMineR::PCA(base_RV,ncp = nb_scores,scale.unit = stand_,graph = FALSE)
  for(j in c(1:nb_scores)){
    score_1<-ACP$ind$coord[,j]
    vect_k<-c(20:L)
    Graphics_estimators_gamma(series = score_1,vect_k = vect_k,
                              Title_graphic =paste0("Paramètre de forme pour le score ",j," de ",nom_variable," standardisee=",stand_),
                              NB_years = NB_year)
  }
  nom_graph<-ifelse(nom_variable=="T(Surcote)",yes="T(S)",
                    no=nom_variable)
  fnct_k<-function(Obs,k){
    LISTE_p<-fonction_analyse_convergence(Obs=Obs,K= k)
    return(LISTE_p)
  }
  vecteur_k<-seq(50,L,by=1)
  # Travail sur la convergence vers un processus l-Pareto. ----------------------------------------------
  
  vecteur_CONVERG<-sapply(X = vecteur_k,FUN = fnct_k,Obs=base_RV)
  MATRICE_moy_coord<-t(cbind.data.frame(vecteur_CONVERG))
  par(mfrow=c(2,2))
  # TITLE_proj<-expression("First moment of "~Theta[M]~"'s projection for "~nom_graph)
  # TITLE_proj<-do.call("substitute", list(TITLE_proj[[1]], list(nom_graph = nom_graph)))                                                                            
  for(j in c(1:ncol(MATRICE_moy_coord))){
    express_h<-expression("Moment with "~h[j])
    express_h<-do.call("substitute", list(express_h[[1]], list(j = j)))
    plot(MATRICE_moy_coord[,j],type="l",ylab=express_h,
         xlab="Exceedances",cex.lab=1.2)
    abline(v=M1,col="red")
    abline(v=M2,col="red")
    # if(j==4){
    #   mtext(outer=TRUE,text =  TITLE_proj,
    #         line = -2)
    # }
  }
  # mtext(outer=TRUE,text = TITLE_proj,
  #       line = -2)
  par(mfrow=c(1,1))
}


#' Graphics_estimators_gamma
#'
#' @param series : vector(float). Variable studied.
#' @param vect_k : vector(int). Vector of the several number of excesses
#' tried. 
#' @param Title_graphic : str. Title of the graphic. 
#' @param NB_years : int. Number of years. 
#'
#' @return Graphic. Evolution of the moment, the Hill 
#' and the ML estimator of the shape parameter. 
#' @export
#'
#' @examples
Graphics_estimators_gamma<-function(series,vect_k,Title_graphic,NB_years=NULL){
  
  NPY<-length(series)/NB_years
  TUNITS<-paste0(NPY,"/year")
  ML_extRemes<-sapply(X =vect_k ,FUN = function_ML_extRemes,
                      data_d=series,typeML="GP",NB_years=NB_years)
  
  Hill_gamma<-sapply(X = vect_k,FUN = function_estimator_hill,series_=series)
  Moments_gamma<-sapply(X = vect_k,FUN=function_estimator_moment,series_=series)
  ML_gamma<-unlist(ML_extRemes[3,])
  data_estimator<-cbind(Moments_gamma,Hill_gamma,
                        ML_gamma)
  Gam<-as.data.frame(melt(t(data_estimator)))
  
  Gam$Var2<-vect_k[Gam$Var2]
  colnames(Gam)<-c("source","number_excesses","gamma_estimed")
  GGothers<-ggplot(data=Gam,aes(x=number_excesses,y =gamma_estimed,col=source,
                                group=interaction(source),
                                linetype=source))+
    geom_line()+
    ylab(expression(gamma))+
    xlab("Number of exceedances")+
    theme(axis.title=element_text(size=15),
          legend.text=element_text(size=10))+
    labs(linetype="Legend",
         caption=paste0("k going from ",min(vect_k)," to ",max(vect_k)))+
    guides(colour="none")+
    ggtitle(Title_graphic)+
    theme(plot.caption = element_text(size = 10))
  #
  print(GGothers)
}

empirical_extremogram<-function(Matrix_couples,inds_select,Tau){
  PIST_empirique<-sapply(Matrix_couples,FUN = pi_s_t,
                         functional_inds=inds_select,vector_u=Tau)
  return(PIST_empirique)
}


# Fonction_quantile_fixe --------------------------------------------------
#' Title
#'
#' @param Obs : dataframe. Processus RV (regular variation). 
#' @param K : int. Nombre d'individus extrêmes
#'
#' @return
#' @export
#'
#' @examples
fonction_analyse_convergence<-function(Obs,K){
  
  NORMES<-apply(Obs,FUN = calcul_norm_L2,MARGIN = 1)
  pas_x<-1/ncol(Obs)
  indice_k<-order(NORMES,decreasing = TRUE)[K]
  Seuil_L2<-NORMES[indice_k]
  Indices<-which(NORMES>Seuil_L2)
  Conservees<-Obs[Indices,]
  NORMES_cons<-NORMES[Indices]
  FORME_d<-t(t(Conservees)%*%diag(NORMES_cons^(-1)))
  vecteur_temps<-c(1:ncol(Obs))/ncol(Obs)
  #Fonction propre. 
  fonction_propre_j<-function(vecteur_temps,j){
    fnct_par_temps<-function(j,t){
      return(sin(2*pi*t*j))
    }
    vecteur_r<-sapply(vecteur_temps,fnct_par_temps,j=j)
    return(vecteur_r)
  }
  LISTE_convergence<-c()
  for(l in c(1:8)){
    fonction_obtenue<-fonction_propre_j(vecteur_temps =vecteur_temps,j=l )

    # approx de Rieman --------------------------------------------------------
    coordonnees<-(FORME_d%*%fonction_obtenue)*(pas_x)
    LISTE_convergence<-c(LISTE_convergence,mean(abs(coordonnees)))
  }
  
  return(LISTE_convergence)
}

Outils_Threshr_choix<-function(Q1,Q2,NT_ths,variable,N_v){
  # Candidates --------------------------------------------------------------
  qu_threshr<- quantile(variable, probs = seq(Q1, Q2,length.out=NT_ths))
  Resultat_threshr<-threshr::ithresh(data=variable,u_vec =qu_threshr,
                                     n_v = N_v)
  return(Resultat_threshr)
}

Outils_POT_graphique<-function(seuil,Q1,Q2,series,dates,titre_variable){
  # Analyse POT seuils pris ------------------------------------------------
  par(mfrow=c(2,2))
  Min_u<-quantile(series,Q1)
  Max_u<-quantile(series,Q2)
  par(cex.lab=1.5,cex.main=1.6)
  POT::tcplot(series,ask = FALSE,which = 1,u.range = c(Min_u,Max_u))
  abline(v=seuil,col="red")
  POT::tcplot(series,ask = FALSE,which = 2,u.range = c(Min_u,Max_u))
  abline(v=seuil,col="red")
  POT::mrlplot(data = series,u.range = c(Min_u,Max_u))
  abline(v=seuil,col="red")
  #"Analysis for ",
  mtext(paste0(titre_variable),outer=TRUE,line = -1,col="red")
  # With time.
  df_time_variable<-cbind.data.frame(1:length(series),
                                     series)
  colnames(df_time_variable)<-c("time","obs")
  
  df_time_variable$time<-dates

  # Attention au saut comme on passe d'hiver en hiver -----------------------
  df_time_variable$time<-as.POSIXct(df_time_variable$time, format="%d/%m/%Y")
  
  # Diplot se base sur des blocs annuels ------------------------------------
  df_time_variable$time<-lubridate::decimal_date(df_time_variable$time)
  
  resultat<-POT::diplot(data = df_time_variable,u.range = c(Min_u,Max_u),nt=200)
  par(cex.lab=1,cex.main=1)
  Indice_dispersion<-abs(resultat$DI-1)
  ind_min<-which.min(Indice_dispersion)[1]
  Threshr_opt<-resultat$thresh[ind_min]
  abline(v=seuil,col="red")
  abline(h=1,col="black")
  par(mfrow=c(1,1))
}


RL_ggplot_boot_opitz<-function(series,seuil,period_years,NPY,titre,
                    nom_variable,plus_simul=FALSE,series_simul=NULL,liste_ML,
                    B,cols_ggplot=NULL,alpha=0.05, 
                    unit_used){

  # Empirique ---------------------------------------------------------------
  xp2 <- ppoints(n = length(series), a = 0)
  sdat <- sort(series) 
  tf_period<-(-1/log(xp2)[sdat > seuil])/NPY
  value_emp<-sdat[sdat>seuil]
  empirique<-cbind.data.frame(tf_period,
                              value_emp)
  colnames(empirique)<-c("transformed.period","sorted.level")
  rate_exceedance<-round(mean(as.numeric(series>seuil)),2)
  m<-period_years*NPY
  Prop<-m*rate_exceedance
  entree<-1-(1/Prop)
  indice_value_twenty<-which(period_years==20)
  print(indice_value_twenty)
  
  # Estim bootstrap ---------------------------------------------------------
  Estim_bootstrap<-Fnct_bootstrap_param(B = B,nb_sample = length(which(series>seuil)),
                       GPD_liste =liste_ML,p =entree,alpha = alpha)
  Base<-data.frame(Estim_bootstrap$niv_retour)
  Base$borne_inf<-ifelse(Base$borne_inf>0,Base$borne_inf,0)
  rvalue_predicted<-Base$estimateur[indice_value_twenty]
  
  series_tronq<-subset(series_simul,series_simul>seuil)
  Niveaux_empiriques_simul<-as.numeric(quantile(x = series_tronq,entree))
  simul_<-cbind.data.frame(period_years,
                           Niveaux_empiriques_simul)
  colnames(simul_)<-c("p_years","simul_niveau")
  
  # GGplot ------------------------------------------------------------------

  GG_RL<-ggplot(data = Base,aes(x=periods_years,y=estimateur))+
    geom_line()+
    geom_line(linetype=0)+
    annotate("point", x = 20, y = rvalue_predicted,colour = "red", 
             size = 4,shape=3)+
    annotate("text", x = 19, y = rvalue_predicted-0.05, 
             label=as.character(round(rvalue_predicted,2)),colour = "red", 
             size = 4)+
    geom_ribbon(mapping = aes(ymin=borne_inf,ymax=borne_sup,col="confidence_band"),alpha=0.15,
                fill="grey", linetype = "dashed")+
    geom_point(data=empirique,aes(x=transformed.period,
                                  y=sorted.level,col="data"))+
    ylab(paste0("return level (",unit_used,")"))+
    xlab("Period P (years)")  
  GG_RL<-GG_RL+
      geom_point(data=simul_,aes(x=p_years,y=simul_niveau,col="simulations"),pch=17)+
      scale_y_continuous(transform = "log")+
      scale_x_continuous(transform="log")+
      labs(colour="Legend")+
     theme(axis.title=element_text(size=15))
      #ggtitle(paste0(titre, " (log-log plot)"))
  if(is.null(cols_ggplot)==FALSE){
    GG_RL<-GG_RL+
      scale_color_manual(values=cols_ggplot)
  }
  print(GG_RL)
}


#' Title
#'
#' @param series : vector.
#' @param seuil : float. 
#' @param period_years : int. 
#' @param NPY : int. 
#' @param titre :str. 
#' @param nom_variable : str. 
#' @param plus_simul : Bool. 
#'
#' @return RL plot with ggplot. 
#' @export
#'
#' @examples
RL_ggplot<-function(series,seuil,period_years,NPY,titre,
                    nom_variable,plus_simul=FALSE,series_simul=NULL,
                    cols_ggplot=NULL,alpha=0.05, 
                    methode_ci="normal",ylim_opt=NULL, 
                    unit_used){
  rate_exceedance<-round(mean(as.numeric(series>seuil)),2)
  series_extreme<-subset(series,series>seuil)
  modele_ev<-fevd( x =series,threshold = seuil,type="GP",
                  time.units = paste0(NPY,"/year"))
  print(modele_ev$results$par)
  scale<-modele_ev$results$par[[1]]
  shape<-modele_ev$results$par[[2]]
  r<-distillery::ci(modele_ev,type = "parameter")
  print("bornes_gamma")
  q_moins<-r[2,1]
  q_plus<-r[2,3]
  test_ad<-goftest::ad.test(x = series_extreme,
                            null = extRemes::"pevd",
                            threshold=seuil,scale=scale,
                            shape=shape,type="GP")
  p_val<-test_ad$p.value
  print(" p value Anderson Darling test")
  print(p_val)
  indice_value_twenty<-which(period_years==20)
  RL_plot<-plot(modele_ev,type = "rl",
       main=paste0("Niveau de retour par EVA pour ",nom_variable),
       rperiods=c(period_years))
  empirique<-RL_plot$empirical

  # Methode delta par defaut -------------------------------------------------
  #############
  
  Modele_vs_emp<-ci(x = modele_ev,alpha = alpha,
                               return.period=period_years,
                                method=methode_ci)
  Base<-cbind.data.frame(Modele_vs_emp[,1],Modele_vs_emp[,2],Modele_vs_emp[,3])
  rvalue_predicted<-Base[indice_value_twenty,2]
  # cas où négatif ----------------------------------------------------------
  
  colnames(Base)<-c("borne_inf","estimateur","borne_sup")
  Base$borne_inf<-ifelse(Base$borne_inf>0,Base$borne_inf,0)
  
  phrase_caption<-paste0("number_extremes=",
  length(series_extreme),", npy_extreme=",round(NPY*rate_exceedance),
  ", threshold=",round(seuil,1),", alpha=",alpha)
  
  if(plus_simul==TRUE){
    m<-period_years*NPY
    Prop<-m*rate_exceedance
    entree<-1-(1/Prop)

    # analyse extremes simul --------------------------------------------------
    series_tronq<-subset(series_simul,series_simul>seuil)
    
    # plot(-1/log(xp2)[sdat > seuil]/(NPY*rate_exceedance), sdat[sdat >
    #                                          seuil])
    simulations<-extRemes::revd(n =length(series_tronq),threshold = seuil,
                                scale = scale, shape = shape,
                                type = "GP")
    df_simul_versus_tronq<-cbind.data.frame(series_tronq, 
                                            simulations)
    colnames(df_simul_versus_tronq)<-c("tronq","simulations_GPD")
    # GG_dens<-ggplot(data=df_simul_versus_tronq,aes(x=tronq))+
    #   geom_density(alpha=0.20,col="blue")+
    #   geom_density(aes(x=simulations_GPD),alpha=0.20,col="yellow")+
    #   ggtitle("Comparaison simul et loi GPD")
    # print(GG_dens)
    Niveaux_empiriques_simul<-as.numeric(quantile(x = series_tronq,entree))
    indices_pris<-sapply(Niveaux_empiriques_simul , function(x,ref){
      return(which.min(abs(x-ref)))
    },ref=series_tronq)
    # enlever valeurs égales.
    #simul_<-cbind.data.frame(-1/log(xp2)[sdat > seuil]/(NPY*rate_exceedance),
    #                         sdat[sdat >
     #                               seuil])
    simul_<-cbind.data.frame(period_years,Niveaux_empiriques_simul)
    colnames(simul_)<-c("p_years","simul_niveau")

    # phrase caption ----------------------------------------------------------
    phrase_caption<-paste0(phrase_caption,", number_simul=",length(series_simul))
    
  }
  GG_RL<-ggplot(data = Base,aes(x=periods_years,y=estimateur))+
    geom_line()+
    geom_line(linetype=0)+
    annotate("point", x = 20, y = rvalue_predicted,colour = "red", 
             size = 4,shape=3)+
    annotate("text", x = 19, y = rvalue_predicted-0.05, 
             label=as.character(round(rvalue_predicted,2)),colour = "red", 
             size = 4)+
    geom_ribbon(mapping = aes(ymin=borne_inf,ymax=borne_sup,col="confidence_band"),alpha=0.15,
              fill="grey", linetype = "dashed")+
    geom_point(data=empirique,aes(x=transformed.period,
                                  y=sorted.level,col="data"))+
    ylab(paste0("return level (",unit_used,")"))+
    
    xlab("Period P (years)")+
    theme(axis.title=element_text(size=15))

  # Use ylim in limits if scale used ------------------------------------------------
  # is not applied otherwise.
  if(plus_simul==TRUE){
    GG_RL<-GG_RL+
      geom_point(data=simul_,aes(x=p_years,y=simul_niveau,col="simulations"),pch=17)
    
  }
  if(is.null(ylim_opt)==FALSE){
    GG_RL<-GG_RL+
        scale_y_continuous(limits = ylim_opt,transform="log")
  }
  else{
    GG_RL<-GG_RL+
        scale_y_continuous(transform="log")
  }
  GG_RL<-GG_RL+
    scale_x_continuous(transform="log")+
    labs(colour="Legend")
    #ggtitle(paste0(titre, " (log-log plot)"))
  #,caption=phrase_caption
  
  if(is.null(cols_ggplot)==FALSE){
    GG_RL<-GG_RL+
      scale_color_manual(values=cols_ggplot)
    
  }
  
  print(GG_RL)
  return(list("plot_RL"=RL_plot,"modele_ev"=modele_ev))
}

fnct_estim_extremo_resample<-function(B,inds_ext,Tau,Matrix_couples,vector_distances){
  N<-nrow(inds_ext)
  Indices_B<-sample(x = c(1:N),size = B,replace = TRUE)
  Estimation_resample<-empirical_extremogram(Matrix_couples = Matrix_couples,
                                             inds_select =inds_ext[Indices_B,],
                                             Tau = Tau)
  df<-cbind.data.frame(Estimation_resample,vector_distances)
  colnames(df)<-c("val_data","delta")
  
  result_delta<-df %>% group_by(delta) %>% summarise(resample=mean(val_data))
  return(result_delta$resample)
}

# Boostrap paramétrique ---------------------------------------------------
Fnct_bootstrap_param<-function(B,nb_sample,GPD_liste,p,alpha){
  sigma_estim<-GPD_liste[["scale"]]
  shape_estim<-GPD_liste[["shape"]]
  mu<-GPD_liste[["threshr"]]
  
  fnct_simul_GPD<-function(n){
    ech<-runif(n = n)
    return(qgp_craft(ech,sigma=sigma_estim,xi =shape_estim)+mu)
  }
  # Simulations GPD ---------------------------------------------------------
  Simulations_GP<-as.matrix(replicate(n=B,expr = fnct_simul_GPD(nb_sample)))
  fnct_estim_1_ech<-function(echantillon,p){
    modele<-fevd(x = echantillon,threshold = mu,type="GP")
    scale_ech<-as.numeric(modele$results$par[[1]])
    shape_ech<-as.numeric(modele$results$par[[2]])
    quantiles_ech<-sapply(p,FUN=function(x){
                  return(qgp_craft(x,sigma = scale_ech,xi = shape_ech)+mu)}
                  )
    return(quantiles_ech)
  }
  fct_estim_params<-function(echantillon,indice_param){
    modele<-fevd(x = echantillon,threshold = mu,type="GP")
    param<-modele$result$par[[indice_param]]
    return(param)
  }
  Niv_retour<-apply(X = Simulations_GP,MARGIN =2,
                    FUN=fnct_estim_1_ech,p=p)
  if(is.null(nrow(Niv_retour))){
    Niv_retour<-matrix(Niv_retour,nrow = length(Niv_retour),
                       ncol = 1)
  }
  Qinf<-apply(X =Niv_retour,MARGIN = 1,FUN=function(x)
    {return(quantile(x,alpha/2))})
  Qsup<-apply(X =Niv_retour,MARGIN = 1,FUN=function(x)
  {return(quantile(x,1-alpha/2))})
  estimateur<-qgp_craft(x = p,xi = shape_estim,sigma = sigma_estim)+mu
  
  # df ----------------------------------------------------------------------
  df<-cbind.data.frame(estimateur,Qinf,Qsup,p)
  colnames(df)<-c("estimateur","borne_inf","borne_sup","p")
  
  # shape -------------------------------------------------------------------
  Estim_scale<-apply(X = Simulations_GP,MARGIN =2,
                    FUN=fct_estim_params,indice_param=1)
  Estim_shape<-apply(X = Simulations_GP,MARGIN =2,
                     FUN=fct_estim_params,indice_param=2)
  return(list("niv_retour"=df,"shape"=Estim_shape,"scale"=Estim_scale))
}

# Utilisation des lois bivariees ------------------------------------------

#' Title
#'
#' @param s : int. Time. 
#' @param t : int. Second time. 
#' @param quantiles_seuil : vector. Set of tail weight.  
#' @param observations : df. Dataframe of observations
#' @param model_bv : str. Bivariate model used. 
#' @param alpha_init : float. Init parameter of the model
#'
#' @return
#' @export
#'
#' @examples
Fonction_ext_bivariee<-function(s,t,quantiles_seuil,observations, 
                                model_bv,alpha_init=NULL){
  series_s<-observations[,s]
  series_t<-observations[,t]
  # chargement quantiles ----------------------------------------------------
  qs<-quantile(series_s,1-quantiles_seuil[s])
  qt<-quantile(series_t,1-quantiles_seuil[t])
  if(is.null(alpha_init)==FALSE){
    modele_bv<-POT::fitbvgpd(cbind(series_s,series_t),
                             c(qs,qt), 
                             model = model_bv,alpha=alpha_init)
  }
  else{
    modele_bv<-POT::fitbvgpd(cbind(series_s,series_t),
                             c(qs,qt), 
                             model = model_bv)
  }
  
  # estimation --------------------------------------------------------------
  return(modele_bv)
}

#' Title
#'
#' @param Number_realisations : number of jumps to do.
#' @param lambda : parameter of the exponential law.
#' @param number_years : integer. Number of years in the original dataframe
#'
#' @return
#' @export
#'
#' @examples
Generator_PP_wake_up<-function(Number_realisations,lambda){
  Times_jump<-rep(NA,Number_realisations)
  for(j in c(1:Number_realisations)){
    # A(x,y)=0 when y>x+1 so only one possibility -----------------------------
    time_jump_j<-rexp(n =1,rate=lambda)
    # Moment of the jth jump --------------------------------------------------
    Times_jump[j]<-time_jump_j
  }
  return(cumsum(Times_jump))
}
fnct_jump_adapt_value<-function(estims_modele,series_values,NYear){
  Temps<-c()
  Longueur<-length(series_values)
  while(length(Temps)!=Longueur){
    liste_lambda<-sapply(series_values,function(x){(1+(estims_modele[[3]])*(x-estims_modele[[1]])/estims_modele[[2]])**(-1/estims_modele[[3]])})
    temps_simul<-sapply(liste_lambda,function(x){return(rexp(1,rate=x))})
    Temps<-c(Temps,min(temps_simul))
    indice<-which(temps_simul==min(temps_simul))
    series_values<-series_values[-indice]
  }
  return(Temps)
}

#' fonction_simul_HTawn
#'
#' @param x : premiere coordonnee.
#' @param y : seconde coordonnee. 
#' @param seuil_x : seuil pour x.  
#' @param seuil_y : seuil pour y. 
#' @param vecteur_x_reg : nouvelles valeurs pour regression (HTawn)
#'
#' @return
#' @export
#'
#' @examples
fonction_simul_HTawn<-function(x,y,seuil_x,seuil_y,vecteur_x_reg){
  
  P_u<-mean(as.numeric(x<seuil_x))
  P_seuil_Y<-mean(as.numeric(y<seuil_y))
  modele_texmex<-texmex::mex(cbind(x,y),which = 1,mth = c(seuil_x,seuil_y))
  noms_liste<-names(modele_texmex$margins$models)
  params_Y0<-modele_texmex$margins$models$y$par
  sigma_Y0<-exp(params_Y0[1])
  gamma_Y0<-params_Y0[2]
  
  params_epsi<-modele_texmex$margins$models$x$par
  sigma_epsi<-exp(params_epsi[1])
  gamma_epsi<-params_epsi[2]
  simul_epsi_realisations<-1-(1-P_u)*texmex::pgpd(q =vecteur_x_reg,sigma = sigma_epsi, 
                                                  xi = gamma_epsi, u = seuil_x, 
                                                  lower.tail = FALSE)
  Laplace_epsi<-modele_texmex$dependence$margins$p2q(simul_epsi_realisations)
  Errors<-modele_texmex$dependence$Z
  Sample_inds<-sample(1:length(Errors),size = length(vecteur_x_reg), 
                      replace = TRUE)
  Errors_sampled<-Errors[Sample_inds,]
  a_estim<-modele_texmex$dependence$coefficients[1]
  b_estim<-modele_texmex$dependence$coefficients[2]
  Pred_Laplace_Y<-a_estim*(Laplace_epsi)+(Laplace_epsi^b_estim)*Errors_sampled
  Pred_Unif_Y<-modele_texmex$dependence$margins$q2p(Pred_Laplace_Y)
  # Deux cas : au-dessus ou en dessous du seuil -----------------------------
  indices_non_ext_YO<-which(Pred_Unif_Y<P_seuil_Y)
  indices_ext_Y0<-which(Pred_Unif_Y>=P_seuil_Y)
  Tirage_non_ext<-quantile(y,
                           Pred_Unif_Y[indices_non_ext_YO])
  poids_queue<-(1-P_seuil_Y)
  Q<-(Pred_Unif_Y[indices_ext_Y0]-P_seuil_Y)/poids_queue
  Tirage_ext<-qgp_craft(x =Q, 
                        sigma = sigma_Y0, 
                        xi = gamma_Y0)+seuil_y
  value_predicted<-rep(NA,length(vecteur_x_reg))
  value_predicted[indices_non_ext_YO]<-Tirage_non_ext
  value_predicted[indices_ext_Y0]<-Tirage_ext
  return(value_predicted)
}

Analyse_threshd_GPD<-function(dates_taken,data,fonction_threshd,n.dens,name,j_show){
  P_valeur_AD_excedent_GPD<-c()
  P_valeur_KS_excedent_GPD<-c()
  
  # Conserver pour plus tard en mémoire. ------------------------------------
  Vect_Theta<-c()
  Vect_gamma<-c()
  Vect_scale<-c()
  Vect_seuil<-c()
  vect_gamma_Moment<-c()
  plot(c(1:37),fonction_threshd,ylim = c(0,0.30),
       main = paste0("Poids mis sur la queue de distribution pour ",name))
  
  par(mfrow=c(3,3))
  for(t in c(1:37)){
    p_ut<-fonction_threshd[t]
    # ML_loi excedents --------------------------------------------------------
    variable_ech_original<-data[,t]
    seuil_t<-quantile(variable_ech_original,1-p_ut)
    name_graph<-ifelse(name=="Surcote","S",name)
    if(t%in%j_show){
      Outils_POT_graphique(series=variable_ech_original,seuil=seuil_t,Q1=0.50,Q2=0.98,
                           dates=dates_taken,
                           titre_variable="")
    }
    
    kernel_dens<-density(x=variable_ech_original,n=n.dens,bw="nrd0")
    indice<-which.max(seuil_t-kernel_dens$x<0)-1
    valeur_fonc_u<-as.numeric(kernel_dens$y[indice])
    sigma_t<-(p_ut/valeur_fonc_u)
    # estimateur_ML -----------------------------------------------------------
    gamma_t<-get.tail.index.estimate(variable_ech_original,1-p_ut,sigma_t)
    
    QGPDpareto<-function(x){
      # si gamma est proche de 0. 
      if(abs(gamma_t)<10^{-4}){
        part1<-(-1)*log(1-x)*sigma_t
        return(part1+seuil_t)
      }
      numerateur<-(1-x)^(-gamma_t)-1
      dename<-gamma_t
      y<-((numerateur/dename)*sigma_t)+seuil_t
      return(y)
    }
    
    queue_distrib<-subset(variable_ech_original,variable_ech_original>seuil_t)
    N_excedents<-length(queue_distrib)
    # estimateur moment -------------------------------------------------------
    gamma_moment_t<-function_estimator_moment(k =N_excedents,
                                              series_ = variable_ech_original)
    vect_gamma_Moment<-c(vect_gamma_Moment,gamma_moment_t)
    
    # Niveaux de retour -------------------------------------------------------
    Dates_converted<-lubridate::decimal_date(as.POSIXct(dates_taken, format="%d/%m/%Y"))
    Dates_converted<-floor(Dates_converted)
    time.rec<-range(Dates_converted)
    NB_annees<-diff(time.rec)
    nber_events_py<-length(variable_ech_original)/NB_annees
    
    # simulationsGPD<-QGPDpareto(ppoints(n = N_excedents))
    # extRemes::qqplot(queue_distrib,simulationsGPD,xlab = "Quantiles empiriques",ylab="Quantiles théoriques")
    # mtext(paste0("Comparaison de ",name_variable," en t=",t," avec une GPD(",round(seuil_t,1),",",round(sigma_t,1),",",round(gamma_t,1),")"))
    Vect_gamma<-c(Vect_gamma,gamma_t)
    Vect_scale<-c(Vect_scale,sigma_t)
    Vect_seuil<-c(Vect_seuil,seuil_t)
    
    # Theta prédit ------------------------------------------------------------
    Theta_t<-ifelse(gamma_t<0,seuil_t+((-1)*sigma_t/gamma_t),Inf)
    Vect_Theta<-c(Vect_Theta,Theta_t)
    
    TEST_AD_pareto_excedents<-goftest::ad.test(queue_distrib,extRemes::"pevd",threshold=seuil_t,scale=sigma_t,shape=gamma_t,type="GP")
    TEST_KS_pareto_excedents<-ks.test(queue_distrib,extRemes::"pevd",threshold=seuil_t,scale=sigma_t,shape=gamma_t,type="GP")
    P_valeur_AD_excedent_GPD<-c(P_valeur_AD_excedent_GPD,TEST_AD_pareto_excedents$p.value)
    P_valeur_KS_excedent_GPD<-c(P_valeur_KS_excedent_GPD,TEST_KS_pareto_excedents$p.value)
  }
  par(mfrow=c(1,1))
  plot(c(1:37),vect_gamma_Moment,type="l",main="Estimateur de gamma (moment) par temps")
  df_test_excedent_GPD<-cbind.data.frame(P_valeur_KS_excedent_GPD,P_valeur_AD_excedent_GPD)
  colnames(df_test_excedent_GPD)<-c("KS","AD")
  GG_TEST_GPD<-ggplot2::ggplot(data = df_test_excedent_GPD,aes(x=1:37,y=KS,col="KS"))+
    geom_point()+
    geom_line(aes(y=KS,col="KS"))+
    geom_point(aes(y=AD,col="AD"))+
    geom_line(aes(y=AD,col="AD"))+
    ylim(c(0,1))+
    ggtitle(paste0("p valeur des tests ",name,"(t)|",name,"(t)>u(t)~GPD(u(t)",",","\u03C3","(t),","\u03B3","(t))"))+
    xlab(label = "t")+
    ylab(label="p valeur")+
    geom_hline(yintercept = 0.05,show.legend = TRUE)+
    scale_color_manual("Test",values = c("red","blue"))+
    labs(caption = paste0("n=",nrow(data),", n.dens=",n.dens))
  print(GG_TEST_GPD)
  
  
  # Garder en memoire EV par temps ------------------------------------------
  df_EV_evol<-cbind.data.frame(Vect_seuil,Vect_scale,Vect_gamma,
                               Vect_Theta,fonction_threshd,P_valeur_AD_excedent_GPD)
  colnames(df_EV_evol)<-c("seuil_t","échelle_t","forme_t",
                          "Theta_t","p_u_t","p_val_ADarling")
  
  # Export de la table ------------------------------------------------------
  write.csv(x=df_EV_evol,file=paste0("Work_RVariations/EVA_",name,".csv"))
  return(P_valeur_AD_excedent_GPD)
}

Analyse_Pareto_per_time<-function(data_Pareto,name,n.dens){
  P_valeur<-c()
  P_valeur_AD<-c()
  par(mfrow=c(3,3))
  for(t in c(1:37)){
    variable_t<-data_Pareto[,t]
    TEST_ks_pareto<-ks.test(x = variable_t,extRemes::"pevd",threshold=1,scale=1,shape=1,type="GP")
    P_valeur<-c(P_valeur,TEST_ks_pareto$p.value)
    TEST_AD_pareto<-goftest::ad.test(x = variable_t,extRemes::"pevd",threshold=1,scale=1,shape=1,type="GP")
    P_valeur_AD<-c(P_valeur_AD,TEST_AD_pareto$p.value)
  }
  par(mfrow=c(1,1))
  df_test_P<-cbind.data.frame(P_valeur,P_valeur_AD)
  colnames(df_test_P)<-c("KS","AD")
  GG_TEST<-ggplot2::ggplot(data = df_test_P,aes(x=1:37,y=KS,col="KS"))+
    geom_point()+
    geom_line(aes(y=KS,col="KS"))+
    geom_point(aes(y=AD,col="AD"))+
    geom_line(aes(y=AD,col="AD"))+
    ylim(c(0,1))+
    ggtitle(paste0("p valeur des tests T(",name,")(t)~Pareto(1)"))+
    xlab(label = "t")+
    ylab(label="p valeur")+
    geom_hline(yintercept = 0.05,show.legend = TRUE)+
    scale_color_manual("Test",values = c("red","blue"))+
    labs(caption = paste0("n=",length(variable_t),", n.dens=",n.dens))
  
  print(GG_TEST)
  
}


# Choice Xo ---------------------------------------------------------------
#' Sample_window
#'
#' @param x_simul : float. l(.) value of the simulated residual. 
#' @param x_window : vector(float). l() values of Epsilon(M). 
#' @param y : vector(float). l() values of X(M-1).
#' @param size_window : int. Number of points in the same window.
#'
#' @return int. Index of the database of the chosen X(M-1). 
#' @export
#'
#' @examples
Sample_window<-function(x_simul,x_window,y,size_window){
  
  series_sort<-sort(x_window)
  index_minimum<-which.min(abs(series_sort-x_simul))
  gap<-(series_sort-x_simul)[index_minimum]
  L<-length(series_sort)-index_minimum
  if(L<(size_window/2)){
    beg<-length(series_sort)-size_window+1
    defs_window<-c(beg:
                      length(series_sort))
  }
  else if(index_minimum<size_window/2){
    end<-size_window
    defs_window<-c(1:end)
  }
  else{
    if(gap>0){
      # the nearest is larger.
      beg<-index_minimum-(size_window/2)
      end<-index_minimum+(size_window/2)-1
    }
    else{
      # the nearest is lower.
      beg<-index_minimum-(size_window/2)+1
      end<-index_minimum+(size_window/2)
    }
    defs_window<-c(beg:end)
  }
  #order gives the initial index used
  L_window<-length(defs_window)
  #ponderate by distance, relative gap ?
  # Wk<-abs(series_sort[defs_window]-x_simul)^(-1)/(x_simul)
  # Wk<-Wk/sum(Wk)
  # summary(Wk)
  #,prob = Wk
  Chosen<-sample(c(1:L_window),size=1)
  index_sample_f<-defs_window[Chosen]
  index_sample<-order(x_window)[index_sample_f]
  return(index_sample)
}
#' Function_AR_p
#'
#' @param j : int. Time index.
#' @param obs_eve : vector(float). The values obtained at each time step for the 
#' previous observation X(M-1).  
#' @param epsilon_t_plus1 : vector(float). Simulated residual.
#' @param model_ar_per_t : df. Summary of the arima model at each time t. 
#'
#' @return vector(float). Simulated extreme time series X(M). 
#' @export
#'
#' @examples
Function_AR_p<-function(j,obs_eve,epsilon_t_plus1,model_ar_per_t){
  
  End<-ncol(model_ar_per_t)-1
  vect_regression<-model_ar_per_t[j,c(1:End)]
  Constant<-c(unlist(model_ar_per_t$Intercept)[j]%*%(1-sum(unlist(vect_regression))))
  X_projection<-as.matrix(obs_eve[j])%*%unlist(vect_regression)+Constant
  return(X_projection+epsilon_t_plus1[j])
}

custom_labels <- function(x) {
  sapply(x, function(val) {
    if (is.na(val)) {
      return(NA)
    } else if (val < 1) {
      sprintf("%.1f", val)
    } else {
      sprintf("%.0f", val)
    }
  })
}


RL_ggplot_cond_ext<-function(series,seuil,period_years,NPY,titre,
                    nom_variable,plus_simul=FALSE,series_simul=NULL,
                    cols_ggplot=NULL,alpha=0.05, 
                    methode_ci="normal",ylim_opt=NULL, 
                    unit_used,Individus_exts){
  rate_exceedance<-round(mean(as.numeric(series>seuil)),2)
  series_extreme<-subset(series,series>seuil)
  modele_ev<-fevd( x =series,threshold = seuil,type="GP",
                   time.units = paste0(NPY,"/year"))
  scale<-modele_ev$results$par[[1]]
  shape<-modele_ev$results$par[[2]]
  r<-distillery::ci(modele_ev,type = "parameter")
  print("bornes_gamma")
  q_moins<-r[2,1]
  q_plus<-r[2,3]
  test_ad<-goftest::ad.test(x = series_extreme,
                            null = extRemes::"pevd",
                            threshold=seuil,scale=scale,
                            shape=shape,type="GP")
  p_val<-test_ad$p.value
  print(" p value Anderson Darling test")
  print(p_val)
  indice_value_twenty<-which(period_years==20)
  pdf(NULL)
  RL_plot<-plot(modele_ev,type = "rl",
                main=paste0("Niveau de retour par EVA pour ",nom_variable),
                rperiods=c(period_years))
  dev.off()
  empirique<-RL_plot$empirical
  # Methode delta par defaut -------------------------------------------------
  #############
  
  Modele_vs_emp<-ci(x = modele_ev,alpha = alpha,
                    return.period=period_years,
                    method=methode_ci)
  Base<-cbind.data.frame(Modele_vs_emp[,1],Modele_vs_emp[,2],Modele_vs_emp[,3])
  rvalue_predicted<-Base[indice_value_twenty,2]
  # cas où négatif ----------------------------------------------------------
  
  colnames(Base)<-c("borne_inf","estimateur","borne_sup")
  Base$borne_inf<-ifelse(Base$borne_inf>0,Base$borne_inf,0)
  
  phrase_caption<-paste0("number_extremes=",
                         length(series_extreme),", npy_extreme=",round(NPY*rate_exceedance),
                         ", threshold=",round(seuil,1),", alpha=",alpha)
  
  if(plus_simul==TRUE){
    m<-period_years*NPY
    Prop<-m*rate_exceedance
    entree<-1-(1/Prop)
    
    # analyse extremes simul --------------------------------------------------
    # with simul
    series_tronq<-subset(series_simul,series_simul>seuil)
    #with data
    #series_data_ext<-series[Individus_exts]
    #series_tronq<-subset(series_data_ext,series_data_ext>seuil)
    # Estimated_quantiles<-sapply(entree,
    #                          fonction_melange_cond,
    #                          series_tronq=series_tronq,
    #                          series=series,seuil=seuil,
    #                          inds_exts=Individus_exts)
    Estimation_probs<-sapply(X = Base$estimateur,
                               cdf_cond_value,
                               series_tronq=series_tronq,
                               series=series,seuil=seuil,
                               inds_exts=Individus_exts)
    # Return_obtained_simul ---------------------------------------------------
    C<-NPY*rate_exceedance
    R_simul<-(C*(1-Estimation_probs))^(-1)
    New_return<-cbind.data.frame(R_simul,Base$estimateur)
    #New_return<-cbind.data.frame(period_years,Estimated_quantiles)
    colnames(New_return)<-c("p_years","simul_niveau")
    
    # phrase caption ----------------------------------------------------------
    phrase_caption<-paste0(phrase_caption,", number_simul=",length(series_simul))
    
  }
  GG_RL<-ggplot(data = Base,aes(x=periods_years,y=estimateur))+
    geom_line()+
    geom_line(linetype=0)+
    annotate("point", x = 20, y = rvalue_predicted,colour = "red", 
             size = 4,shape=3)+
    annotate("text", x = 19, y = rvalue_predicted-0.05, 
             label=as.character(round(rvalue_predicted,2)),colour = "red", 
             size = 4)+
    geom_ribbon(mapping = aes(ymin=borne_inf,ymax=borne_sup,col="confidence_band"),alpha=0.15,
                fill="grey", linetype = "dashed")+
    geom_point(data=empirique,aes(x=transformed.period,
                                  y=sorted.level,col="data"))+
    ylab(paste0("return level (",unit_used,")"))+
    
    xlab("Period P (years)")
  
  # Use ylim in limits if scale used ------------------------------------------------
  # is not applied otherwise.
  if(plus_simul==TRUE){
    GG_RL<-GG_RL+
      geom_point(data=New_return,aes(x=p_years,y=simul_niveau,col="simulations"),pch=17)
    
  }
  if(is.null(ylim_opt)==FALSE){
    GG_RL<-GG_RL+
      scale_y_continuous(limits = ylim_opt,transform="log",
                         labels = label_number(accuracy = 0.001))
  }
  else{
    GG_RL<-GG_RL+
      scale_y_continuous(transform="log",
                         labels = label_number(accuracy = 0.001))
  }
  GG_RL<-GG_RL+
    scale_x_continuous(labels = custom_labels,transform = "log")+
    labs(colour="Legend")+
    theme(axis.title=element_text(size=15),
          legend.text=element_text(size=10))
   ## ggtitle(paste0(titre, " (log-log plot)"))
  
  if(is.null(cols_ggplot)==FALSE){
    GG_RL<-GG_RL+
      scale_color_manual(values=cols_ggplot)
  }
  
  print(GG_RL)
  return(list("plot_RL"=RL_plot,"modele_ev"=modele_ev))
}
cdf_cond_value<-function(quantile_emp,series_tronq,series,seuil,inds_exts){
  # inds_exts = extreme time series (L2 norm) --------------------------------
  denom<-length(which(series>seuil))
  # ext ---------------------------------------------------------------------
  prop_ext_estim<-mean(as.numeric((series_tronq<quantile_emp)))
  
  # intersect extreme TS AND above threshold --------------------------------
  
  prop_a<-length(intersect(inds_exts,which(series>seuil)))/denom
  moitie<-prop_ext_estim*prop_a
  
  # non_ext -----------------------------------------------------------------
  inds_non_exts<-c(1:length(series))[-inds_exts]
  
  # intersect non extreme TS AND above threshold ----------------------------
  
  case_2<-length(intersect(which((series>seuil)&(series<quantile_emp)),
                           inds_non_exts))/denom
  return(moitie+case_2)
}

cdf_melange_cond<-function(quantile_emp,series_tronq,series,seuil,inds_exts,u_given){
  
  # inds_exts = extreme time series (L2 norm) --------------------------------
  denom<-length(which(series>seuil))
  # ext ---------------------------------------------------------------------
  prop_ext_estim<-mean(as.numeric((series_tronq<quantile_emp)))

  # intersect extreme TS AND above threshold --------------------------------

  prop_a<-length(intersect(inds_exts,which(series>seuil)))/denom
  moitie<-prop_ext_estim*prop_a
  
  # non_ext -----------------------------------------------------------------
  inds_non_exts<-c(1:length(series))[-inds_exts]

  # intersect non extreme TS AND above threshold ----------------------------

  case_2<-length(intersect(which((series>seuil)&(series<quantile_emp)),
                           inds_non_exts))/denom
  return(moitie+case_2-u_given)
}

fonction_melange_cond<-function(u_given,series_tronq,series,seuil,inds_exts){

  # find a coherent level of quantile level ---------------------------------

  # knowing what we give as entry (u_given) -------------------------------------------

  
  return(uniroot(f = cdf_melange_cond,interval = c(0,10),
                 series_tronq=series_tronq,
                 u_given=u_given,inds_exts=inds_exts,
                 seuil=seuil,series=series
                 )$root)
}


#' Resamples_tendencies
#'
#' @param extremes_indus : observed extreme time series. 
#' @param B : int. Number of resamples. 
#' @param list_Q : vector(float). Vector of quantiles analysed
#'
#' @return vector(float). Percent value found at each time
#' @export
#'
#' @examples
Resamples_tendencies<-function(extremes_indus,B,list_Q){
  N<-nrow(extremes_indus)
  Indexes_B<-sample(x = c(1:N),size = B,replace = TRUE)
  ech<-extremes_indus[Indexes_B,]
  LIST_results<-matrix(NA,ncol=ncol(ech),nrow = length(list_Q))
  for(l in c(1:length(list_Q))){
    Q<-list_Q[l]
    result_q<-apply(ech,function(x){return(quantile(x,Q))},
                               MARGIN=2)
    LIST_results[l,]<-result_q
  }
  return(LIST_results)
}


# Functions_graphic -------------------------------------------------------

# Short-cut for graphics using the scores/variance ------------------------

function_expression_prop_variance_j<-function(j,prop_variance){
  first_object<-expression(C[j])
  first_object<-do.call("substitute", 
                        list(first_object[[1]], list(j = j)))
  expr_j<-expression(first_object~" ("~prop1~"% of the variance)")
  return_obj<-do.call("substitute", 
                      list(expr_j[[1]], list(prop1 = prop_variance[j],
                                             first_object=first_object)))
  return(return_obj)
}
