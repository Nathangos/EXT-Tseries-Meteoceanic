#' Generateur_beta
#'
#' @param u : float. Realisation of an uniform law. 
#' @param quant : float. Threshold in the uniform scale.
#' @param params_lois : list. Params of the normal law followed for the bulk.  
#' @param seuil_ev : float. Threshold of the extreme law.
#' @param ech_ev : float. Scale parameter of the extreme law. 
#' @param gam_ev : float. Shape parameter of the extreme law. 
#'
#' @return Float. Realisation of the law of beta.
#' @export
#'
#' @examples
Generateur_beta<-function(u,quant,params_lois,seuil_ev,ech_ev,gam_ev){
  
  if(u<quant){
    return(qnorm(p = u,mean = params_lois[[1]],
                 sd=params_lois[[2]]))
  }
  else{
    qprime<-(u-quant)/(1-quant)
    return(extRemes::qevd(p=qprime,
                          threshold =seuil_ev,shape = gam_ev,
                          scale=ech_ev,type="GP"))
  }
  
}


# Generation functions ----------------------------------------------------
#' Generation_curves_poly
#'
#' @param beta_ : float. Beta parameter. 
#' @param c 
#' @param L_courbe 
#' @param begg 
#' @param end 
#' @param dephase 
#'
#' @return
#' @export
#'
#' @examples
Generation_curves_poly<-function(beta_,c,L_courbe,begg,end,dephase){
  
  B<-1
  vect_time<-seq.int(from = begg,to = end,length.out = L_courbe)
  Mid_<-(begg+end)/2
  sig<-(end-begg)/6
  if(dephase==TRUE){
    delta<-rnorm(n = 1,sd = sig-1/10)
    Mid_<-Mid_+delta
  }
  alpha_b<-2*B*Mid_
  Poly<--B*(vect_time)^2+alpha_b*vect_time+beta_
  Curve<-Poly
  Sin_part<-sin(pi/2*(seq.int(from = begg,to = end,length.out = L_courbe))*beta_)+Mid_
  return(Curve+Sin_part)
}


#' Generation_curves_normal_LAW
#'
#' @param beta_ : float. Inverse of the standard deviation of the normal density.
#' @param dephase : Bool. Indicate if the mean changes randomly.
#' @param begg : float. Minimum value of the time variable.
#' @param end : float. Maximum value of the time variable. 
#' @param L_courbe : float. Number of measures for one time series. 
#'
#' @return Evolution of the normal density knowing the value of beta.
#' @export
#'
#' @examples
Generation_curves_normal_LAW<-function(beta_,dephase,begg,end,L_courbe){
  
  delta<-0
  if(dephase==TRUE){
    #6 au depart
    delta<-(rnorm(1)/L_courbe)*4
  }
  # change the moment of the mean.
  mu<-(end+begg)/2+delta
  vect_time<-seq.int(from = begg,to = end,length.out = L_courbe)
  # abs and sign are used to have at the end centered distribution at each time.
  return(beta_+sign(beta_)*(abs(beta_)/sqrt(2*pi))*exp(-((beta_)^2/2)*(vect_time-mu)^2))
  #return(beta*vect_time**2)
}

Generation_curves_sin<-function(beta_,dephase,begg,end,L_courbe){
  delta<-0
  if(dephase==TRUE){
    delta<-(rnorm(1,sd=2)/L_courbe)*pi
  }
  mid<-(pi/2)+delta
  seq_temps<-pi/2*(seq.int(from = begg,to = end,length.out = L_courbe))/beta_+mid
  return(sapply(X =seq_temps,FUN = function(x){return(sin(x))}))
}

#' f_marginales_Pareto
#'
#' @param variable : vector. Realisations of the law at time t
#' @param p_u : float. Weight put on the distribution tail.
#' @param n.dens : int. Parameter of the kernel density estimator.
#'
#' @return list. EV parameters and density estimated, observation 
#' converted in the Pareto scale
#' @export
#'
#' @examples
f_marginales_Pareto<-function(variable,p_u,n.dens){
  kernel_dens<-density(x=variable,n=n.dens)
  threshold<-as.numeric(quantile(variable,1-p_u))
  extremes<-subset(variable,variable>threshold)
  index<-which.max(threshold-kernel_dens$x<0)-1
  value_fonc_u<-as.numeric(kernel_dens$y[index])
  scale_fonc<-(p_u/value_fonc_u)
  k<-which.max(variable-threshold<0)-1
  gamma<-get.tail.index.estimate(variable,1-p_u,scale_fonc)
  gap<-diff(kernel_dens$x)[1]
  # Riemann
  F_dens<-cumsum(kernel_dens$y)*gap
  non_ext<-subset(variable,variable<=threshold)
  intervals<-c(-Inf,kernel_dens$x+gap/2)
  interval_vecteur<-findInterval(non_ext,vec = intervals)
  n<-length(variable)
  
  indices_Excedent<-which(variable>threshold)
  vector_unif<-rep(NA,length(variable))
  Excesses<-variable[indices_Excedent]-threshold
  
  Fcraft<-function(x,xi,sigma){
    return(pgp_craft(x=x,xi=xi,sigma=sigma))
  }
  Unif_excesses<-(1-p_u)+(p_u)*sapply(Excesses,Fcraft,
                                      sigma = scale_fonc,xi = gamma)
  vector_unif[-indices_Excedent]<-F_dens[interval_vecteur]
  vector_unif[indices_Excedent]<-Unif_excesses
  
  
  vector_pareto<-(1-vector_unif)^(-1)
  return(list("obs"=as.vector(vector_pareto),"kernel_dens"=kernel_dens,"gamma"=gamma,
              "scale"=scale_fonc,"threshold"=threshold,"p_u"=p_u))
  
}

#' f_marginales_all_Pareto
#'
#' @param indice : int. Time index. 
#' @param data_to_tf : dataframe. Temporal series to convert. 
#' @param p_u : float. Weight put on the extremes at each time. 
#' @param n.dens : int. Parameter of the kernel density estimator. 
#'
#' @return
#' @export
#'
#' @examples
f_marginales_all_Pareto<-function(indice,data_to_tf,p_u,n.dens){
  
  variable<-data_to_tf[,indice]
  p_UI<-p_u[indice]
  return(f_marginales_Pareto(variable = variable,p_u = p_UI,n.dens=n.dens))
}

# Functions reconversion --------------------------------------------------
#' functions_reconversion_Pareto
#'
#' @param K : list(list). For each given t, the estimated density of the time series 
#' at time t. 
#' @param variable_uplift : vector. One simulated individual at the Pareto scale. 
#' @param list_evt : list. Parameters of the extreme law for each value of t. 
#'
#' @return
#' @export
#'
#' @examples
function_reconversion_Pareto<-function(K,variable_uplift,list_evt){
  L_dim<-length(variable_uplift)
  individu_ext<-sapply(1:L_dim,function_reconv_each_dim_Pareto,K=K,v_diff=variable_uplift,
                       list_evt=list_evt)
  return(individu_ext)
}
#' fonction_reconv_each_dim_Pareto
#'
#' @param index_dimension : int. Time index. 
#' @param K : list(list). For each given t, the estimated density of the time series 
#' at time t. 
#' @param v_diff : vector. Simulated vector at the Pareto scale. 
#' @param list_evt : list. Parameters of the extreme law for each value of t.
#'
#' @return Value in the original scale at time index_dimension for the individual 
#'         simulated. 
#' @export
#'
#' @examples
function_reconv_each_dim_Pareto<-function(index_dimension,K,v_diff,list_evt){
  
  valeur_x<-v_diff[index_dimension]
  kernel_DENS<-K[[index_dimension]]
  p_u<-list_evt[["p_u"]]
  p_Ui<-p_u[[index_dimension]]
  scale_fonc<-list_evt[["scale"]][[index_dimension]]
  gamma<-list_evt[["gamma"]][[index_dimension]]
  u_f<-list_evt[["threshold"]][[index_dimension]]
  indice_S<-which.max(u_f-kernel_DENS$x<0)-1
  step<-diff(kernel_DENS$x)[1]
  integrale<-cumsum(kernel_DENS$y)*step
  variable_unif<-(-(valeur_x)^(-1))
  var_reference<-variable_unif+1
  if(-p_Ui>=variable_unif){
    indice_d_min<-which.max(var_reference-integrale<0)-1
    valeur<-as.numeric(kernel_DENS$x[indice_d_min])
    return(valeur)
  }
  else{
    value_p<-(variable_unif+p_Ui)/p_Ui
    evalue<-qgp_craft(x = value_p,sigma = scale_fonc,xi = gamma)+u_f
    return(as.numeric(evalue))
  }
}

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

##### Calculate the L2 norm of a time series. #####

#' calcul_norme_L2
#'
#' @param series 
#'
#' @return int. L2 norm of the series. 
#' @export
#'
#' @examples
calcul_norme_L2<-function(series){
  L<-length(series)
  pas<-pas<-(1/(L-1))
  
  # Indicatrice de positivité -----------------------------------------------
  Sigma<-sum(subset(series,series>0)**2)
  norme_carree<-Sigma*pas
  return(norme_carree**(1/2))
}

#' function_reconstitution_trajectory_std
#'
#' @param Vector_coords : dataframe. Pca coordinates obtained.
#' @param Base_functions_p : dataframe. Eigen vectors obtained. 
#' @param NB_dim : int. Chosen number of PCA vectors. 
#' @param mu_t : vector(float). Mean function. 
#' @param var_t : vector(float). Standard value function.
#' @param BOOL_STD : boolean.  Type of series used (standardised or not).
#'
#' @return Angular components of the observations. 
#' @export
#'
#' @examples
#' 
#' 
function_reconstitution_trajectory_std<-function(Vector_coords,Base_functions_p,NB_dim,mu_t,sd_t,BOOL_STD){
  L<-nrow(Vector_coords)
  Matrix_score<-diag(Vector_coords[,1],ncol = L)
  feigen1<-Base_functions_p[,1]
  Matrice_repmfonction<-matrix(rep(feigen1,L),ncol=L)
  Simulations<-t(Matrice_repmfonction%*%Matrix_score)
  Mat_mean<-matrix(rep(mu_t,L),byrow = TRUE,nrow=L)
  if(NB_dim==1){
    if(BOOL_STD==TRUE){
      Shape_ACP<-Simulations%*%diag(sd_t)+Mat_mean
    }
    else{
      Shape_ACP<-Simulations+Mat_mean
    }
  }
  else{
    for(j in c(2:NB_dim)){
      MS<-diag(Vector_coords[,j] ,ncol=L)
      Matrice_f_p2<-matrix(rep(Base_functions_p[,j],L),ncol=L)
      Simulations<-Simulations+t(Matrice_f_p2%*%MS)
    }
    if(BOOL_STD==TRUE){
      Shape_ACP<-Simulations%*%diag(sd_t)+Mat_mean
    }
    else{
      Shape_ACP<-Simulations+Mat_mean
    }
  }
  
  return(Shape_ACP)
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

fnct_select_colonne<-function(df,numero){
  return(df[numero,])
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
#' function_Structure_Matrice
#'
#' @param NB_dim : int. Number of coordinates
#'
#' @return Generate the matrix used in Rvine for the structure
#' @export
#'
#' @examples
function_Structure_Matrice<-function(NB_dim){
  
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
    theme(axis.title=element_text(size=15))+
    labs(linetype="Legend",
         caption=paste0("k going from ",min(vect_k)," to ",max(vect_k)))+
    guides(colour="none")+
    ggtitle(Title_graphic)
  #
  print(GGothers)
}
### We use here Frechet margins so we will check that we are above 1 at each time. 
Function_margin_check<-function(series){
  resultat_test<-all(series>0)
  return(resultat_test)
}

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

#' fonction_ML_extRemes
#'
#' @param k 
#' @param donnees 
#' @param typeML
#'
#' @return extRemes package estimator of the shape parameter with its boundaries. 
#' @export
#'
#' @examples
fonction_ML_extRemes<-function(k,donnees,typeML="GP",NB_annees){
  
  ## Ne fonctionne pas si le gamma est inférieur à -1/2.
  # Méthode paramétrique.
  donnees_ord<-sort(donnees,decreasing = TRUE)
  seuil<-donnees_ord[k]
  L<-length(donnees)/NB_annees
  theta_k<-extRemes::fevd(donnees,threshold = seuil,type=typeML,
                          time.units = paste0(L,"/year"))
  
  tryCatch(expr={
    r<-distillery::ci(theta_k,type = "parameter")
    # formule est 
    # std_<-sqrt(inverse_hessienne)[2,2]
    # ga-qnorm(0.975)*std_
    # ga+qnorm(0.975)*std_
    q_moins<-r[2,1]
    estimateur<-r[2,2]
    q_plus<-r[2,3]
    liste_r<-list("q_moins"=q_moins,"q_plus"=q_plus,"estimateur"=estimateur)
    return(liste_r)
  },
  error=function(e){
    estimateur<-theta_k$results$par[[2]]
    liste_r<-list("q_moins"=NA,"q_plus"=NA,"estimateur"=estimateur)
    return(liste_r)}
  )
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


#' function_analyse_convergence
#'
#' @param Obs : dataframe. Processus RV (regular variation). 
#' @param K : int. Nombre d'individus extrêmes
#'
#' @return
#' @export
#'
#' @examples
function_analyse_convergence<-function(Obs,K){
  
  NORMES<-apply(Obs,FUN = calcul_norme_L2,MARGIN = 1)
  pas_x<-1/ncol(Obs)
  index_k<-order(NORMES,decreasing = TRUE)[K]
  Threshold_L2<-NORMES[index_k]
  Indices<-which(NORMES>Threshold_L2)
  Conserved<-Obs[Indices,]
  NORMES_cons<-NORMES[Indices]
  Shape_d<-t(t(Conserved)%*%diag(NORMES_cons^(-1)))
  vect_time<-c(1:ncol(Obs))/ncol(Obs)
  #Eigen function 
  Eigen_function_j<-function(vect_time,j){
    fnct_per_time<-function(j,t){
      return(sin(2*pi*t*j))
    }
    vect_r<-sapply(vect_time,fnct_per_time,j=j)
    return(vect_r)
  }
  LIST_convergence<-c()
  for(l in c(1:8)){
    function_obtained<-Eigen_function_j(vect_time =vect_time,j=l )
    
    # approx de Rieman --------------------------------------------------------
    Coordinates<-(Shape_d%*%function_obtained)*(pas_x)
    LIST_convergence<-c(LIST_convergence,mean(abs(Coordinates)))
  }
  
  return(LIST_convergence)
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
  index_sample_f<-defs_window[sample(c(1:L_window),size=1)]
  index_sample<-order(x_window)[index_sample_f]
  return(index_sample)
}

#' Function_AR_p
#'
#' @param j : int. Time index.
#' @param obs_eve : vector(float). The values obtained at each time step for the 
#' previous observation X(M-1).  
#' @param epsilon_t_plus1 : vector(float). Simulated residual.
#' @param Modele_ar_par_t : df. Summary of the arima model at each time t. 
#'
#' @return vector(float). Simulated extreme time series X(M). 
#' @export
#'
#' @examples
Function_AR_p<-function(j,obs_eve,epsilon_t_plus1,Modele_ar_par_t){
  
  End<-ncol(Modele_ar_par_t)-1
  vect_regression<-Modele_ar_par_t[j,c(1:End)]
  Constant<-c(unlist(Modele_ar_par_t$Intercept)[j]%*%(1-sum(unlist(vect_regression))))
  X_projection<-obs_eve[j]%*%unlist(vect_regression)+Constant
  return(X_projection+epsilon_t_plus1[j])
}



# ML--fonctions -----------------------------------------------------------

#' fnct_ML_realite_simulations
#'

#' @param ech_train : list. Training elements with simulated time series (or values)
#'                    and recorded time series (or values).
#' @param ech_test : list. Testing elements with simulated time series (or values)
#'                    and recorded time series (or values).
#' @param type_noyau : str. Chosen kernel.
#' @return
#' @export
#'
#' @examples
fnct_ML_realite_simulations<-function(ech_train,ech_test,type_noyau){

  Modif<-FALSE
  simul_train<-ech_train$simul
  simul_test<-ech_test$simul
  real_train<-ech_train$real
  real_test<-ech_test$real

  if(is.null(nrow(simul_train))){
    simul_train<-matrix(simul_train,byrow = FALSE,ncol = 1,
                        nrow = length(simul_train))
    simul_test<-matrix(simul_test,byrow = FALSE,ncol = 1,
                       nrow = length(simul_test))
    real_train<-matrix(real_train,byrow = FALSE,ncol = 1,
                       nrow = length(real_train))
    real_test<-matrix(real_test,byrow = FALSE,ncol = 1,
                      nrow = length(real_test))
  }
  realisations<-rbind(real_train,simul_train)
  mean_realisations<-apply(X = realisations,MARGIN = 2,
                           FUN = mean)
  sd_realisations<-apply(X = realisations,MARGIN = 2,
                         FUN = sd)
  realisations<-as.data.frame(scale(as.matrix(realisations)))
  # Variable cible : indicatrice simulation ---------------------------------

  # Echantillon d'apprentissage ---------------------------------------------
  # y=1 for simulations et y=0 for observations
  realisations$y<-c(rep(0,nrow(real_train)),rep(1,nrow(simul_train)))
  realisations$y<-as.factor(realisations$y)
  realisations<-realisations[sample(nrow(realisations)),]
  modele_svm<-e1071::svm(y~.,data=realisations,kernel=type_noyau,
                         probability=TRUE)

  predictions<-modele_svm$fitted
  DF_test<-rbind(real_test,simul_test)
  DF_test<-as.data.frame(scale(DF_test,center=mean_realisations,
                               scale=sd_realisations))
  rownames(DF_test)<-c(1:nrow(DF_test))
  y<-c(rep(0,nrow(real_test)),rep(1,nrow(simul_test)))
  # predictions--SVM ---------------------------------------------------------
  predictions_modele<-predict(modele_svm,DF_test,
                              decision.values = TRUE)
  y_pred<-as.numeric(as.character(predictions_modele))
  scores<-attr(predictions_modele,"decision.values")
  probs_fonction_chosen<-as.numeric((1+exp((-1)*scores))^(-1))

  # Take into account the order ---------------------------------------------
  # sign(score) gives yi.  --------------------------------------------------
  if(y_pred[1]!=round(probs_fonction_chosen[1])){
    probs_fonction_chosen<-(1-probs_fonction_chosen)
  }
  return(list("Modèle"=modele_svm,"resultat_conf_train"=table(realisations$y,predictions),
              "resultat_conf_test"=table(y,predictions_modele),
              "type_noyau"=type_noyau,
              "predictions_modele"=predictions_modele,
              "probs_y"=probs_fonction_chosen,
              "y"=y_pred,
              "y_test"=y))
}


# RandomForest ------------------------------------------------------------

#' fnct_RF
#'
#' @param ech_train : list. Training elements with simulated time series (or values)
#'                    and recorded time series (or values).
#' @param ech_test : list. Testing elements with simulated time series (or values)
#'                    and recorded time series (or values).
#' @param NB_trees : int. Number of trees used in the forest.
#'
#' @return
#' @export
#'
#' @examples
fnct_RF<-function(ech_train,ech_test,NB_trees){

  Modif<-FALSE
  simul_train<-ech_train$simul
  simul_test<-ech_test$simul
  real_train<-ech_train$real
  real_test<-ech_test$real

  if(is.null(nrow(simul_train))){
    simul_train<-matrix(simul_train,byrow = FALSE,ncol = 1,
                        nrow = length(simul_train))
    simul_test<-matrix(simul_test,byrow = FALSE,ncol = 1,
                       nrow = length(simul_test))
    real_train<-matrix(real_train,byrow = FALSE,ncol = 1,
                       nrow = length(real_train))
    real_test<-matrix(real_test,byrow = FALSE,ncol = 1,
                      nrow = length(real_test))
  }
  realisations<-rbind(real_train,simul_train)
  mean_realisations<-apply(X = realisations,MARGIN = 2,
                           FUN = mean)
  sd_realisations<-apply(X = realisations,MARGIN = 2,
                         FUN = sd)
  realisations<-as.data.frame(scale(as.matrix(realisations)))
  # Variable cible : indicatrice simulation ---------------------------------

  # Echantillon d'apprentissage ---------------------------------------------
  # y=1 for simulations et y=0 for observations
  L<-ncol(realisations)
  realisations$y<-c(rep("0",nrow(real_train)),rep("1",nrow(simul_train)))
  realisations$y<-as.factor(realisations$y)
  realisations<-realisations[sample(nrow(realisations)),]
  Input_rf<-realisations[,1:L]
  if(is.null(nrow(Input_rf))){
    Input_rf<-matrix(Input_rf,byrow = FALSE,ncol = 1,
                     nrow = length(Input_rf))
  }
  modele_rf<-randomForest::randomForest(Input_rf,y=realisations$y,
                                        data=realisations,
                                        ntree=NB_trees)
  predictions<-modele_rf$predicted
  DF_test<-rbind(real_test,simul_test)
  if(ncol(DF_test)==1){
    DF_test<-scale(DF_test,center=mean_realisations,
                   scale=sd_realisations)
  }
  else{
    DF_test<-as.data.frame(scale(DF_test,center=mean_realisations,
                                 scale=sd_realisations))
  }
  y<-c(rep(0,nrow(real_test)),rep(1,nrow(simul_test)))
  # predictions--RF ---------------------------------------------------------
  probs_fonction_chosen<-predict(modele_rf,DF_test,
                                 type = "prob")[,2]
  y_pred<-round(probs_fonction_chosen)
  return(list("Modèle"=modele_rf,"resultat_conf_train"=table(realisations$y,predictions),
              "resultat_conf_test"=table(y,y_pred),
              "probs_y"=probs_fonction_chosen,
              "y"=y_pred,
              "y_test"=y))
}

#' fnct_logistic_glm
#'
#' @param ech_train : list. Training elements with simulated time series (or values)
#'                    and recorded time series (or values).
#' @param ech_test : list. Testing elements with simulated time series (or values)
#'                    and recorded time series (or values).
#' @param fnct_lien : str. Link function used. 
#'
#' @return List. Fitted model with main results. 
#' @export
#'
#' @examples
fnct_logistic_glm<-function(ech_train,ech_test,fnct_lien){
  
  Modif<-FALSE
  simul_train<-ech_train$simul
  simul_test<-ech_test$simul
  real_train<-ech_train$real
  real_test<-ech_test$real
  
  if(is.null(nrow(simul_train))){
    simul_train<-matrix(simul_train,byrow = FALSE,ncol = 1,
                        nrow = length(simul_train))
    simul_test<-matrix(simul_test,byrow = FALSE,ncol = 1,
                       nrow = length(simul_test))
    real_train<-matrix(real_train,byrow = FALSE,ncol = 1,
                       nrow = length(real_train))
    real_test<-matrix(real_test,byrow = FALSE,ncol = 1,
                      nrow = length(real_test))
  }
  realisations<-rbind(real_train,simul_train)
  mean_realisations<-apply(X = realisations,MARGIN = 2, 
                           FUN = mean)
  sd_realisations<-apply(X = realisations,MARGIN = 2, 
                         FUN = sd)
  realisations<-as.data.frame(scale(as.matrix(realisations)))
  # Variable cible : indicatrice simulation ---------------------------------
  
  # Echantillon d'apprentissage ---------------------------------------------
  
  # y=1 pour les simulations et y=0 pour les observations
  realisations$y<-c(rep(0,nrow(real_train)),rep(1,nrow(simul_train)))
  realisations$y<-as.factor(realisations$y)
  
  # Shuffle -----------------------------------------------------------------
  realisations<-realisations[sample(nrow(realisations)),]
  modele_logit<-glm(y~.,data=realisations,family = binomial(link = fnct_lien))
  predictions<-round(modele_logit$fitted)
  
  # Echantillon_test --------------------------------------------------------
  y<-c(rep(0,nrow(real_test)),rep(1,nrow(simul_test)))
  vect<-c(real_test,simul_test)
  DF_test<-rbind(real_test,simul_test)
  DF_test<-as.data.frame(scale(DF_test,center=mean_realisations, 
                               scale=sd_realisations))
  rownames(DF_test)<-c(1:nrow(DF_test))
  # predictions--logit ---------------------------------------------------------
  probs_model<-as.numeric(predict(modele_logit,DF_test,
                                  type="response"))
  predictions_modele<-round(probs_model)
  return(list("Modèle"=modele_logit,"resultat_conf_train"=table(realisations$y,predictions),
              "resultat_conf_test"=table(y,predictions_modele),
              "type_fnct_lien"=fnct_lien,"y_test"=y,"probs_y"=probs_model))
}

#' Function_extract_elt
#'
#' @param modele_stat : model. Return object of craft GLM, SVM and RF functions.
#'
#' @return Key performance values.
#' @export
#'
#' @examples
Function_extract_elt<-function(modele_stat){
  
  # ROC ---------------------------------------------------------------------
  reality_<-modele_stat$y_test
  probs_algorithm_<-modele_stat$probs_y
  curve_roc_<-pROC::roc(reality_,
                        probs_algorithm_,plot=FALSE, 
                        quiet=TRUE)
  index_mid_<-which.min(abs(curve_roc_$thresholds-0.5))
  Liste_ROC_<-list("curve"=curve_roc_, 
                   "index_mid"=index_mid_)
  Mat_confTest_<-modele_stat$resultat_conf_test
  df_<-matrix(0,nrow=2,ncol=2)
  if(ncol(Mat_confTest_)!=2){
    indice<-as.numeric(colnames(Mat_confTest_))+1
    df_[,indice]<-Mat_confTest_
  }
  else{
    df_<-matrix(c(Mat_confTest_),nrow = 2,ncol=2)
  }
  Num<-2*(det(df_))
  Sum_lignes<-prod(apply(df_,MARGIN = 1,FUN = function(x){return(sum(x))}))
  
  Sum_cols<-prod(apply(df_,MARGIN = 2,FUN = function(x){return(sum(x))}))
  Denom<-Sum_lignes+Sum_cols
  
  # HSS near 0= no skill ----------------------------------------------------
  
  HSS_score<-Num/Denom
  Prop_<-t(apply(df_,MARGIN = 1,
                 function(x){return(x/sum(x))}))
  R_<-sum(diag(df_))/sum(df_)
  
  
  return(list("1"=Prop_,"2"=R_,"3"=Liste_ROC_,"4"=HSS_score))
}

#' fnct_calcul_correct_proportion
#'
#' @param Base_simul : dataframe of simulated time series (or L2 norm).
#' @param Base_data : dataframe of observed time series (or L2 norm).
#' @param hyp_param : vector(str). Link function and kernel used.
#' @param K : int. K fold parameter.
#'
#' @return
#' @export
#'
#' @examples
fnct_calcul_correct_proportion<-function(Base_simul, Base_data,hyp_param,
                                         type_sampling,K){
  if(!is.null(nrow(Base_data))){
    Nb_ech<-nrow(Base_data)
    Nb_simul<-nrow(Base_simul)
  }
  else{
    Nb_ech<-length(Base_data)
    Nb_simul<-length(Base_simul)
  }
  Base_ind_s<-c(1:Nb_simul)
  Base_ind_r<-c(1:Nb_ech)

  # sampling of real observations ------------------------------------------------------
  Indices_train_r<-sample(Base_ind_r,size = round(Nb_ech*0.7),replace = FALSE)
  Indices_test_r<-which(!(Base_ind_r%in% Indices_train_r))

  # Sampling of simulations, size  ------------------------------------------
  # because size simul >> size data.
  Indices_s<-sample(Base_ind_s,size =Nb_ech,replace = FALSE)
  if(!is.null(nrow(Base_simul))){
    Sub_Simul<-Base_simul[Indices_s,]
  }
  else{
    Sub_Simul<-Base_simul[Indices_s]
  }
  if(type_sampling=="simple"){
    Inds_s_test<-sample(x = c(1:Nb_ech),size = round(Nb_ech*0.3))
    Inds_r_test<-sample(x = c(1:Nb_ech),size = round(Nb_ech*0.3))
    # Construction of the df --------------------------------------------------
    if(!is.null(nrow(Base_data))){
      ech_train<-list("simul"=Sub_Simul[-Inds_s_test,],
                      "real"=Base_data[-Inds_r_test,])
      ech_test<-list("simul"=Sub_Simul[Inds_s_test,],
                     "real"= Base_data[Inds_r_test,])
    }
    else{
      ech_train<-list("simul"=Sub_Simul[-Inds_s_test],
                      "real"=Base_data[-Inds_r_test])
      ech_test<-list("simul"=Sub_Simul[Inds_s_test],
                     "real"= Base_data[Inds_r_test])
    }

    modele_SVM<-fnct_ML_realite_simulations(ech_train = ech_train,ech_test=ech_test,
                                            type_noyau = hyp_param[1])
    SVM_results<-Function_extract_elt(modele_stat = modele_SVM)
    Prop_SVM<-SVM_results[[1]]
    R_SVM<-SVM_results[[2]]
    Liste_ROC_SVM<-SVM_results[[3]]
    HSS_SVM<-SVM_results[[4]]

    modele_GLM<-fnct_logistic_glm(fnct_lien = hyp_param[2],
                                    ech_test=ech_test,ech_train=ech_train)
    GLM_results<-Function_extract_elt(modele_stat = modele_GLM)
    Prop_GLM<-GLM_results[[1]]
    R_GLM<-GLM_results[[2]]
    Liste_ROC_GLM<-GLM_results[[3]]
    HSS_<-GLM_results[[4]]



    # RF ----------------------------------------------------------------------
    modele_RF<-fnct_RF(ech_train = ech_train,ech_test=ech_test,
                       NB_trees = as.numeric(hyp_param[3]))
    RF_results<-Function_extract_elt(modele_stat = modele_RF)
    Prop_RF<-RF_results[[1]]
    R_RF<-RF_results[[2]]
    Liste_ROC_RF<-RF_results[[3]]
    HSS_rf<-RF_results[[4]]

  }

  return(list("GLM"=list("accuracy"=R_GLM,
                         "modeles"=modele_GLM,
                         "ROC"=Liste_ROC_GLM,
                         "HSS"=HSS_),
              "SVM"=list("accuracy"=R_SVM,
                         "modeles"=modele_SVM,
                         "ROC"=Liste_ROC_SVM,
                         "HSS"=HSS_SVM),
              "RF"=list("accuracy"=R_RF,
                        "modeles"=modele_RF,
                        "ROC"=Liste_ROC_RF,
                        "HSS"=HSS_rf)))
}

#' Running_perfs_ML
#'
#' @param Base_simul : dataframe (or vector). Simulated time series (or values).
#' @param Base_data : dataframe (or vector). Observed time series (or values).
#' @param hyp_param : list. Hyper parameters such as the link function (GLM), the kernel (SVM)
#'                    and the number of trees (RF).
#' @param NB_times : int. Number of repetitions of the experiment (=number of models used).
#' @param alpha_prop : float. Confidence level used.
#' @param K : int. Not used. Number of blocks used in cross validation.
#' @param type_sampling : str. Type of experiment used (simple or CV).
#' @param title_ROC : str. Titre used in ROC graphics.
#' @param NB_shown_ROC : str. Number of ROC curves shown.
#'
#' @return
#' @export
#'
#' @examples
Running_perfs_ML<-function(Base_simul, Base_data,hyp_param,NB_times,alpha_prop,K,
                           type_sampling, title_ROC,NB_shown_ROC){


  Resultat<-replicate(n = NB_times,fnct_calcul_correct_proportion(Base_simul=Base_simul,
                                                                  Base_data = Base_data,
                                                                  hyp_param =hyp_param,type_sampling=type_sampling,
                                                                  K=K))
  Echs<-sample(c(1:NB_times),size = NB_shown_ROC)
  Whole_curve<-apply(Resultat,MARGIN = 2,
                     function(x){return(x$SVM$ROC$curve)})

  GGROC_several<-pROC::ggroc(Whole_curve[Echs],alpha=0.4)+
    guides(col="none")+
    geom_abline(slope = 1, intercept = 1, linetype = "dashed", color = "black")+
    ggtitle(paste0("SVM (", hyp_param[1],") ROC curves for ",title_ROC," (shown=",NB_shown_ROC,", calculated=",NB_times,")"))
  print(GGROC_several)
  Whole_curve_GLM<-apply(Resultat,MARGIN = 2,
                         function(x){return(x$GLM$ROC$curve)})
  GGROC_several_2<-pROC::ggroc(Whole_curve_GLM[Echs],alpha=0.4)+
    geom_abline(slope = 1, intercept = 1, linetype = "dashed", color = "black")+
    guides(col="none")+
    ggtitle(paste0("GLM (", hyp_param[2],") ROC curves for ",title_ROC," (shown=",NB_shown_ROC,", calculated=",NB_times,")"))
  print(GGROC_several_2)

  Whole_curve_RF<-apply(Resultat,MARGIN = 2,function(x){return(x$RF$ROC$curve)})
  GGROC_several_RF<-pROC::ggroc(Whole_curve_RF[Echs],alpha=0.4)+
    geom_abline(slope = 1, intercept = 1, linetype = "dashed", color = "black")+
    guides(col="none")+
    ggtitle(paste0("Random Forest ROC curves for ",title_ROC," (shown=",NB_shown_ROC,", calculated=",NB_times,")"))
  print(GGROC_several_RF)

  # Accuracy ----------------------------------------------------------------

  GLM_ACCURACY<-apply(Resultat,MARGIN = 2,function(x){return(x$GLM$accuracy)})
  SVM_ACCURACY<-apply(Resultat,MARGIN = 2,function(x){return(x$SVM$accuracy)})
  RF_ACCURACY<-apply(Resultat,MARGIN = 2,function(x){return(x$RF$accuracy)})
  Accuracy<-cbind(SVM_ACCURACY,GLM_ACCURACY,RF_ACCURACY)
  Q1<-apply(X = Accuracy,MARGIN = 2,FUN = function(x){return(round(quantile(x,(alpha_prop/2)),2)*100)})
  Q2<-apply(X = Accuracy,MARGIN = 2,FUN = function(x){return(round(quantile(x,1-(alpha_prop/2)),2)*100)})
  df<-cbind.data.frame(Q1,Q2)
  resultat<-cbind.data.frame(c("SVM","GLM","RF"),paste0(df$Q1,"-",df$Q2))
  colnames(resultat)<-c("Modele","qu_accuracy")

  # HSS ----------------------------------------------------------------
  GLM_HSS<-apply(Resultat,MARGIN = 2,function(x){return(x$GLM$HSS)})
  SVM_HSS<-apply(Resultat,MARGIN = 2,function(x){return(x$SVM$HSS)})
  RF_HSS<-apply(Resultat,MARGIN = 2,function(x){return(x$RF$HSS)})
  HSS<-cbind(SVM_HSS,GLM_HSS,RF_HSS)
  Q1_h<-apply(X = HSS,MARGIN = 2,FUN = function(x){return(round(quantile(x,(alpha_prop/2)),2)*100)})
  Q2_h<-apply(X = HSS,MARGIN = 2,FUN = function(x){return(round(quantile(x,1-(alpha_prop/2)),2)*100)})
  df<-cbind.data.frame(Q1_h,Q2_h)
  resultat_h<-cbind.data.frame(c("SVM","GLM","RF"),paste0("(",df$Q1,")-(",df$Q2,")"))
  colnames(resultat_h)<-c("Modele","qu_hss")

  return(list("Accuracy"=resultat,
              "HSS"=resultat_h,
              "Ensemble_r"=Resultat))
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
  indices_B<-sample(x = c(1:N),size = B,replace = TRUE)
  ech<-extremes_indus[indices_B,]
  LIST_results<-matrix(NA,ncol=ncol(ech),nrow = length(list_Q))
  for(l in c(1:length(list_Q))){
    Q<-list_Q[l]
    result_q<-Percent_Q<-apply(ech,function(x){return(quantile(x,Q))},
                                        MARGIN=2)
    LIST_results[l,]<-result_q
  }
  return(LIST_results)
}

empirical_extremogram<-function(Matrix_couples,inds_select,Tau){
  PIST_empirique<-sapply(Matrix_couples,FUN = pi_s_t,
                         functional_inds=inds_select,vector_u=Tau)
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

fnct_estim_extremo_resample<-function(B,inds_ext,Tau,Matrix_couples,vector_distances){
  N<-nrow(inds_ext)
  Indices_B<-sample(x = c(1:N),size = B,replace = TRUE)
  Estimation_resample<-empirical_extremogram(Matrix_couples = Matrice_couples,
                                          inds_select =inds_ext[Indices_B,],
                                          Tau = Tau)
  df<-cbind.data.frame(Estimation_resample,vecteur_distances)
  colnames(df)<-c("val_data","delta")
  
  result_delta<-df %>% group_by(delta) %>% summarise(resample=mean(val_data))
  return(result_delta$resample)
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
Analyse_extreme_proj<-function(base_RV,nb_scores,stand_,name_variable,L,NB_year,M1,M2){
  
  ACP<-FactoMineR::PCA(base_RV,ncp = nb_scores,scale.unit = stand_,graph = FALSE)
  for(j in c(1:nb_scores)){
    score_1<-ACP$ind$coord[,j]
    vect_k<-c(20:L)
    Graphics_estimators_gamma(series = score_1,vect_k = vect_k,
                              Title_graphic =paste0("Shape parameter for the score ",j," of ",name_variable," standardised=",stand_),
                              NB_years = NB_year)
  }
  name_graph<-ifelse(name_variable=="T(Surcote)",yes="T(S)",
                    no=name_variable)
  fnct_k<-function(Obs,k){
    LISTE_p<-function_analyse_convergence(Obs=Obs,K= k)
    return(LISTE_p)
  }
  vect_k<-seq(50,L,by=1)

  vect_CONVERG<-sapply(X = vect_k,FUN = fnct_k,Obs=base_RV)
  MATRIX_moy_coord<-t(cbind.data.frame(vect_CONVERG))
  par(mfrow=c(2,2))
  TITLE_proj<-expression("First moment of "~Theta[M]~"'s projection for "~name_graph)
  TITLE_proj<-do.call("substitute", list(TITLE_proj[[1]], list(name_graph = name_graph)))                                                                            
  for(j in c(1:ncol(MATRIX_moy_coord))){
    express_h<-expression("Moment with "~h[j])
    express_h<-do.call("substitute", list(express_h[[1]], list(j = j)))
    plot(MATRIX_moy_coord[,j],type="l",ylab=express_h,
         xlab="Exceedances")
    abline(v=M1,col="red")
    abline(v=M2,col="red")
    if(j==4){
      mtext(outer=TRUE,text =  TITLE_proj,
            line = -2)
    }
  }
  mtext(outer=TRUE,text = TITLE_proj,
        line = -2)
  par(mfrow=c(1,1))
}

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
  q_moins<-r[2,1]
  q_plus<-r[2,3]
  test_ad<-goftest::ad.test(x = series_extreme,
                            null = extRemes::"pevd",
                            threshold=seuil,scale=scale,
                            shape=shape,type="GP")
  p_val<-test_ad$p.value
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
    series_tronq<-subset(series_simul,series_simul>seuil)
    Estimated_quantiles<-sapply(entree,
                                fonction_melange_cond,
                                series_tronq=series_tronq,
                                series=series,seuil=seuil,
                                inds_exts=Individus_exts)
    
    # Return_obtained_simul ---------------------------------------------------
    #C<-NPY*rate_exceedance
    #R_simul<-(C*(1-Estimation_probs))^(-1)
    New_return<-cbind.data.frame(period_years,Estimated_quantiles)
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
      geom_point(data=New_return,aes(x=p_years,y=simul_niveau,col="simul"),pch=17)
    
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
    theme(axis.title=element_text(size=15))
  ## ggtitle(paste0(titre, " (log-log plot)"))
  
  if(is.null(cols_ggplot)==FALSE){
    GG_RL<-GG_RL+
      scale_color_manual(values=cols_ggplot)
  }
  
  print(GG_RL)
  return(list("plot_RL"=RL_plot,"modele_ev"=modele_ev))
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
  
  return(uniroot(f = cdf_melange_cond,interval = c(0,max(series_tronq)),
                 series_tronq=series_tronq,
                 u_given=u_given,inds_exts=inds_exts,
                 seuil=seuil,series=series
  )$root)
}