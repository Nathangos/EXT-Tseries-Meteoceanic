source("fonctions/fonctions.R")

# Prendre en compte seuil changeant ---------------------------------------
f_marginales_uniformeV2<-function(variable,p_u_dim,n.dens){
  kernel_dens<-density(x=variable,n=n.dens)
  seuil<-as.numeric(quantile(variable,1-p_u_dim))
  indice<-which.max(seuil-kernel_dens$x<0)-1
  valeur_fonc_u<-as.numeric(kernel_dens$y[indice])
  scale_fonc<-(p_u_dim/valeur_fonc_u)
  gamma<-get.tail.index.estimate(variable,1-p_u_dim,scale_fonc)
  theta<-extRemes::fevd(variable,threshold = seuil,type = "GP",method = "MLE")
  scale_foncP<-theta$results$par[[1]]
  gammaP<-theta$results$par[[2]]
  print("#######################")
  print("# Nouvelle dimension --------------------------------------#")
  print("Selon le package")
  print(paste("Forme vaut",gammaP,", Scale=",scale_foncP))
  print("Selon les fonctions de l'article")
  print(paste("Forme vaut",gamma,", Scale=",scale_fonc))
  print("#######################")

  # indices_moins
  ecart<-diff(kernel_dens$x)[1]
  # Riemann
  F_dens<-cumsum(kernel_dens$y)*ecart
  non_ext<-subset(variable,variable<=seuil)
  intervals<-c(-Inf,kernel_dens$x+ecart/2)
  interval_vecteur<-findInterval(non_ext,vec = intervals)
  interval_vecteur[length(interval_vecteur)] = Inf
  # nouvelles valeurs
  n<-length(variable)
  indices_Excedent<-which(variable>seuil)
  vecteur_unif<-rep(NA,length(variable))
  vecteur_unif[variable<=seuil]<-F_dens[interval_vecteur]
  #vecteur_unif[variable>seuil]<-(1-p_u_dim)+(p_u_dim)*extRemes::pevd(variable[indices_Excedent],threshold =seuil,scale=scale_fonc,shape=gamma,type="GP")
  # Methode de github. 
  vecteur_unif[indices_Excedent]<-(1-p_u_dim)+(p_u_dim)*pgp_craft(variable[indices_Excedent]-seuil,sigma = scale_fonc,xi = gamma)
  vecteur_unif_decalage<-vecteur_unif-1
  return(list("obs"=as.vector(vecteur_unif_decalage),"kernel_dens"=kernel_dens))
}
f_marginales_uniforme_EVOLTAU<-function(indice,donnees,p_u,n.dens){
  variable<-donnees[,indice]
  p_UI<-p_u[indice]
  return(f_marginales_uniformeV2(variable = variable,p_u_dim = p_UI,n.dens=n.dens))
}
#' fonction_uplifting_EVOLTAU
#'
#' @param indice : int. Indice dans la base de l'observation.
#' @param donnees_transf : df. Donnees obtenues apres conversion uniforme. 
#' @param bornes_unif : vector. Bornes de la variable uniforme.
#' @param p_u : vector. Vecteur des poids mis dans les queues de distribution.
#' @param K : List. Liste compose des objets density en chaque temps.
#' @param donnees_originales : df. Donnees brutes. 
#'
#' @return Vector. L'individu extreme transforme en l'indice . 
#' @export
#'
#' @examples
fonction_uplifting_EVOLTAU<-function(indice,donnees_transf,vecteur_V,p_u,K,donnees_originales,type_hissage){
  variable_transformee<-donnees_transf[indice,]
  Norme<-fonc_norm_inv(variable_transformee)
  Vprime<-vecteur_V[indice]
  s<-Vprime/Norme
  Indices_dim<-1:length(variable_transformee)
  X_uplift<-sapply(Indices_dim,fonction_uplift_dim,s=s,p_u=p_u,K=K,xext=variable_transformee,donnees_originales=donnees_originales,type_hissage)
  return(X_uplift)
}

#' fonction_uplift_dim
#'
#' @param indice_dim : int. Dimension traitee. 
#' @param s : float. Nouvelle echelle a donner. 
#' @param p_u : vector. Vecteur des poids mis dans les queues de distribution.
#' @param K : List. Liste compose des objets density en chaque temps.
#' @param xext : Individu extreme considere.
#' @param donnees_originales : df. Observations brutes. 
#'
#' @return Float. Valeur en la dimensiona pres elevation/abaissement. 
#' @export
#'
#' @examples
fonction_uplift_dim<-function(indice_dim,s,p_u,K,xext,donnees_originales,type_hissage){
  x_obs_dim<-xext[indice_dim]
  # Selection du bon quantile. 
  ref_densite<-donnees_originales[,indice_dim]
  kernel_dens<-K[[indice_dim]]
  p_Ui<-p_u[indice_dim]
  seuil<-quantile(ref_densite,1-p_Ui)
  pas<-diff(kernel_dens$x)[1]
  F_x<-cumsum(kernel_dens$y)*pas
  rang_ds_noyau<-which.max(seuil-kernel_dens$x<0)-1
  nom_h<-names(type_hissage)[1]
  if(nom_h=="seuils"){
    u_marg<-F_x[rang_ds_noyau]-1
  }
  else{
    u_marg<-type_hissage[[nom_h]]
  }
  #
  # Se demander pour chaque marginale si l'observation est d'un côté ou de l'autre de la distribution. 
  # Uplift premier cas.
  if(x_obs_dim>=u_marg){
    x_uplift<-as.numeric(s*x_obs_dim)
  }
  else{
    rapport<-((1+s*u_marg)/(1+u_marg))
    x_uplift<-as.numeric(-1+ rapport*(1+x_obs_dim))
  }
  return(x_uplift)
}
#' fonction_reconversion_EVOLTAU
#'
#' @param K : List. Liste compose des objets density en chaque temps.
#' @param variable_uplift : Vector. Sortie de fonction_uplifting_EVOLTAU.
#' @param p_u : vector. Vecteur des poids mis dans les queues de distribution.
#' @param donnees_originales : df. Observations brutes.  
#'
#' @return Observation extreme reconverti.
#' @export
#'
#' @examples
fonction_reconversion_EVOLTAU<-function(K,variable_uplift,p_u,donnees_originales){
  L_dim<-length(variable_uplift)
  INDU_ext_post_reconversion<-sapply(1:L_dim,fonction_reconversion_par_dim_EVOLTAU,K=K,v_diff=variable_uplift,p_u=p_u,donnees_originales=donnees_originales)
  return(INDU_ext_post_reconversion)
}
#' fonction_reconversion_par_dim_EVOLTAU
#'
#' @param indice_dimension : int. Indice de la dimension.  
#' @param K : List. Liste compose des objets density en chaque temps.
#' @param v_diff : Vector. Observation a reconvertir.
#' @param p_u : vector. Vecteur des poids mis dans les queues de distribution.
#' @param donnees_originales : df. Observations brutes.
#'
#' @return Float. Valeur de l'obervation extreme reconvertie en indice_dimension. 
#' @export
#'
#' @examples
fonction_reconversion_par_dim_EVOLTAU<-function(indice_dimension,K,v_diff,p_u,donnees_originales){
  valeur_x<-v_diff[indice_dimension]
  kernel_DENS<-K[[indice_dimension]]
  p_Ui<-p_u[indice_dimension]
  colonne<-donnees_originales[,indice_dimension]
  u_f<-as.numeric(quantile(colonne,1-p_Ui))
  indice_S<-which.max(u_f-kernel_DENS$x<0)-1
  step<-diff(kernel_DENS$x)[1]
  integrale<-cumsum(kernel_DENS$y)*step
  var_lag<-valeur_x+1
  if(-p_Ui>=valeur_x){
    # Si observation non extreme, on utilise l'estimateur empirique.
    indice_d_min<-which.max(var_lag-integrale<0)-1
    valeur<-as.numeric(kernel_DENS$x[indice_d_min])
    return(valeur)
  }
  else{
    # Si on a une observation extreme, on doit utiliser la deuxieme fonction 
    # de repartition. 
    valeur_fonc_u<-as.numeric(kernel_DENS$y[indice_S])
    scale_fonc<-(p_Ui/valeur_fonc_u)
    gamma<-get.tail.index.estimate(colonne,1-p_Ui,scale_fonc)
    valeur_p<-(valeur_x+p_Ui)/p_Ui
    if(valeur_p>1){
      print(indice_dimension)
    }
    #evalue<-extRemes::qevd(p = valeur_p,loc = u_f,scale = scale_fonc,shape = gamma,type="GP")
    # Le deuxieme ne prend pas en compte l'ajout du seuil.
    evalue<-qgp_craft(x = valeur_p,sigma = scale_fonc,xi = gamma)+u_f
    return(as.numeric(evalue))
  }
}

