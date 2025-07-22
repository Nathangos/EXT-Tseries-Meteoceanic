#source("../fonctions/fonctions.R")

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



# Fonctions reconversion --------------------------------------------------
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


pi_s_t<-function(individus_fonctionnel,vecteur_u,couples_indices){
  t<-couples_indices[1]
  s<-couples_indices[2]
  numerateur<-sum(as.numeric((individus_fonctionnel[,s]>vecteur_u[s])&(individus_fonctionnel[,t]>vecteur_u[t])))
  denominateur<-sum(as.numeric(individus_fonctionnel[,s]>vecteur_u[s]))
  return(numerateur/denominateur)
}
variogramme_BESSEL<-function(t,s,kappa,tau){
  h<-t-s
  h_norme<-abs(h/tau)
  return(kappa*(1-h_norme))
}
sigma_BESSEL<-function(s,kappa,tau){
  h_norme<-abs(s/tau)
  return(kappa*(1-h_norme))
}

